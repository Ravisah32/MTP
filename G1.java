import java.io.*;
import java.util.*;

public class G1 {

    public static void main(String[] args) throws IOException {
        if (args.length < 5) {
            System.out.println(
                    "Usage: java BFSFairRange_MultiType_Fast <csvFile> <startMin> <startMax> <numTypes> <epsilon>");
        } else {
            String csvFile = args[0];
            double startMin = Double.parseDouble(args[1]);
            double startMax = Double.parseDouble(args[2]);
            int numTypes = Integer.parseInt(args[3]);
            int epsilon = Integer.parseInt(args[4]);

            List<Person> people = readCSV(csvFile);

            // sort people by salary
            people.sort(Comparator.comparingDouble(p -> p.salary));
            double[] salaries = new double[people.size()];
            for (int i = 0; i < people.size(); i++)
                salaries[i] = people.get(i).salary;

            // map categories to index
            Set<String> categorySet = new HashSet<>();
            for (Person p : people)
                categorySet.add(p.category);
            List<String> categories = new ArrayList<>(categorySet);
            Map<String, Integer> catIndex = new HashMap<>();
            for (int i = 0; i < categories.size(); i++)
                catIndex.put(categories.get(i), i);

            // build prefix counts per category
            int C = categories.size();
            int N = people.size();
            int[][] prefix = new int[C][N + 1];
            for (int i = 0; i < N; i++) {
                int c = catIndex.get(people.get(i).category);
                for (int j = 0; j < C; j++) {
                    prefix[j][i + 1] = prefix[j][i] + (j == c ? 1 : 0);
                }
            }

            Range startRange = new Range(startMin, startMax);

            System.out.println("Counts in original range:");
            printCategoryCounts(salaries, prefix, categories, startRange);

            long queryStart = System.currentTimeMillis();
            Range fairRange = bfsFindFirstFairRange(salaries, prefix, categories, startRange, numTypes, epsilon);
            long queryEnd = System.currentTimeMillis();

            System.out.println("Query executed in " + (queryEnd - queryStart) + " ms");

            if (fairRange != null) {
                System.out.println("Fair range found: [" + fairRange.min + ", " + fairRange.max + "]");
                System.out.println("Counts in final fair range:");
                printCategoryCounts(salaries, prefix, categories, fairRange);
            } else {
                System.out.println("No fair range found.");
            }
        }
    }

    // ======== Data structures ========

    static class Person {
        double salary;
        String category;

        Person(double s, String c) {
            this.salary = s;
            this.category = c;
        }
    }

    static class Range {
        double min, max;

        Range(double min, double max) {
            this.min = min;
            this.max = max;
        }
    }

    // PQ node: carries bounds by index for precision + similarity to original
    static class Node {
        int L, R; // inclusive indices into sorted salaries
        double sim; // Jaccard vs original range output

        Node(int L, int R, double sim) {
            this.L = L;
            this.R = R;
            this.sim = sim;
        }
    }

    // ======== IO helpers ========

    static List<Person> readCSV(String filename) throws IOException {
        List<Person> people = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
            String line = br.readLine(); // header
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty())
                    continue;
                // Expect: salary,category
                String[] parts = splitCSV(line);
                if (parts.length < 2)
                    continue;
                double sal = Double.parseDouble(parts[0].trim());
                String cat = parts[1].trim();
                people.add(new Person(sal, cat));
            }
        }
        return people;
    }

    // minimal CSV split (handles simple quoted commas)
    static String[] splitCSV(String line) {
        List<String> out = new ArrayList<>();
        StringBuilder sb = new StringBuilder();
        boolean inQ = false;
        for (int i = 0; i < line.length(); i++) {
            char ch = line.charAt(i);
            if (ch == '"')
                inQ = !inQ;
            else if (ch == ',' && !inQ) {
                out.add(sb.toString());
                sb.setLength(0);
            } else
                sb.append(ch);
        }
        out.add(sb.toString());
        return out.toArray(new String[0]);
    }

    static void printCategoryCounts(double[] salaries, int[][] prefix, List<String> categories, Range r) {
        int L = lowerBound(salaries, r.min);
        int R = upperBound(salaries, r.max) - 1; // inclusive
        if (L > R || L < 0 || R >= salaries.length) {
            System.out.println("No salaries in range.");
            return;
        }
        for (int c = 0; c < categories.size(); c++) {
            int cnt = prefix[c][R + 1] - prefix[c][L];
            System.out.println(categories.get(c) + ": " + cnt);
        }
    }

    // ======== Core BFS (PriorityQueue) ========

    static Range bfsFindFirstFairRange(double[] salaries,
            int[][] prefix,
            List<String> categories,
            Range startRange,
            int numTypes,
            int epsilon) {
        int n = salaries.length;
        if (n == 0)
            return null;

        // translate user's startRange to index range [L0..R0]
        int L0 = lowerBound(salaries, startRange.min);
        int R0 = upperBound(salaries, startRange.max) - 1;
        if (L0 < 0)
            L0 = 0;
        if (R0 >= n)
            R0 = n - 1;
        if (L0 > R0)
            return null; // empty original range

        // original output set cardinality (indices form the set elements)
        final int origSize = (R0 - L0 + 1);

        // Precompute a helper to get counts quickly (for fairness and Jaccard)
        // We will compute |A∩B| from index intervals directly.

        // PriorityQueue: max-heap by similarity; tiebreak by smaller width, then
        // lexicographically
        PriorityQueue<Node> pq = new PriorityQueue<>((a, b) -> {
            int cmp = Double.compare(b.sim, a.sim);
            if (cmp != 0)
                return cmp;
            int wa = a.R - a.L, wb = b.R - b.L;
            if (wa != wb)
                return Integer.compare(wa, wb);
            if (a.L != b.L)
                return Integer.compare(a.L, b.L);
            return Integer.compare(a.R, b.R);
        });

        // visited ranges by index bounds
        boolean[][] visitedRowCompressed = null;
        // fallback: use a hash set of keys "L#R"
        Set<Long> visited = new HashSet<>();

        // push the original node with sim=1.0
        pq.add(new Node(L0, R0, 1.0));
        visited.add(packKey(L0, R0));

        // If original is already fair, return immediately
        if (isFair(salaries, prefix, L0, R0, categories.size(), numTypes, epsilon)) {
            return new Range(salaries[L0], salaries[R0]);
        }

        while (!pq.isEmpty()) {
            Node cur = pq.poll();

            // Check fairness for this node
            if (isFair(salaries, prefix, cur.L, cur.R, categories.size(), numTypes, epsilon)) {
                return new Range(salaries[cur.L], salaries[cur.R]);
            }

            // Generate neighbors: expand/shrink left/right to next distinct boundary
            // Each neighbor should be a minimal move to alter membership.
            for (int[] nb : neighborsByIndex(salaries, cur.L, cur.R)) {
                int nl = nb[0], nr = nb[1];
                if (nl < 0 || nr >= n || nl > nr)
                    continue;
                long key = packKey(nl, nr);
                if (visited.contains(key))
                    continue;
                visited.add(key);

                // Compute Jaccard similarity w.r.t. original [L0..R0]
                // Intersection size between intervals = overlap length if any
                int inter = intervalIntersectionSize(L0, R0, nl, nr);
                int candSize = nr - nl + 1;
                int union = origSize + candSize - inter;
                double sim = (union == 0) ? 1.0 : ((double) inter) / union;

                pq.add(new Node(nl, nr, sim));
            }
        }

        return null; // not found
    }

    // === Fairness check (plug your original logic here if different) ===
    // Uses per-category counts on [L..R]. A simple parity-style example:
    // - At least 'numTypes' categories present in the range
    // - And max difference between any two category counts <= epsilon
    static boolean isFair(double[] salaries, int[][] prefix, int L, int R,
            int C, int numTypes, int epsilon) {
        int present = 0;
        int minCnt = Integer.MAX_VALUE, maxCnt = Integer.MIN_VALUE;
        for (int c = 0; c < C; c++) {
            int cnt = prefix[c][R + 1] - prefix[c][L];
            if (cnt > 0)
                present++;
            minCnt = Math.min(minCnt, cnt);
            maxCnt = Math.max(maxCnt, cnt);
        }
        if (present < numTypes)
            return false;
        // If all zero (empty range) => not fair
        if (maxCnt == 0)
            return false;
        return (maxCnt - minCnt) <= epsilon;
    }

    // ======== Neighbor generation on index bounds ========
    // Returns a list of [L', R'] neighbors by minimally moving a boundary
    static List<int[]> neighborsByIndex(double[] s, int L, int R) {
        int n = s.length;
        List<int[]> v = new ArrayList<>(4);

        // Move left boundary right (shrink-left) to next distinct salary
        int nL = nextRightDistinctLeft(s, L, R);
        if (nL != -1)
            v.add(new int[] { nL, R });

        // Move left boundary left (expand-left) to previous distinct salary
        int pL = prevLeftDistinctLeft(s, L);
        if (pL != -1)
            v.add(new int[] { pL, R });

        // Move right boundary left (shrink-right) to previous distinct salary
        int nR = prevLeftDistinctRight(s, L, R);
        if (nR != -1)
            v.add(new int[] { L, nR });

        // Move right boundary right (expand-right) to next distinct salary
        int pR = nextRightDistinctRight(s, R);
        if (pR != -1)
            v.add(new int[] { L, pR });

        return v;
    }

    // Helpers to step to next/prev distinct boundary indices

    // Move L → right to the first index > current group of equal salaries
    static int nextRightDistinctLeft(double[] s, int L, int R) {
        int i = L;
        while (i <= R && s[i] == s[L])
            i++;
        return (i <= R) ? i : -1;
        // If equal block consumes to R, shrinking-left would empty -> avoid
    }

    // Move L → left to include the previous block of equal salaries
    static int prevLeftDistinctLeft(double[] s, int L) {
        if (L == 0)
            return -1;
        int i = L - 1;
        double val = s[i];
        while (i - 1 >= 0 && s[i - 1] == val)
            i--;
        return i;
    }

    // Move R → left to drop the current block of equal salaries
    static int prevLeftDistinctRight(double[] s, int L, int R) {
        int j = R;
        while (j >= L && s[j] == s[R])
            j--;
        return (j >= L) ? j : -1;
    }

    // Move R → right to include the next block of equal salaries
    static int nextRightDistinctRight(double[] s, int R) {
        int n = s.length;
        if (R == n - 1)
            return -1;
        int j = R + 1;
        double val = s[j];
        while (j + 1 < n && s[j + 1] == val)
            j++;
        return j;
    }

    // ======== Index / search helpers ========

    // first index i with s[i] >= x
    static int lowerBound(double[] s, double x) {
        int l = 0, r = s.length;
        while (l < r) {
            int m = (l + r) >>> 1;
            if (s[m] < x)
                l = m + 1;
            else
                r = m;
        }
        return l;
    }

    // first index i with s[i] > x
    static int upperBound(double[] s, double x) {
        int l = 0, r = s.length;
        while (l < r) {
            int m = (l + r) >>> 1;
            if (s[m] <= x)
                l = m + 1;
            else
                r = m;
        }
        return l;
    }

    static int intervalIntersectionSize(int aL, int aR, int bL, int bR) {
        int L = Math.max(aL, bL);
        int R = Math.min(aR, bR);
        return (L <= R) ? (R - L + 1) : 0;
    }

    static long packKey(int L, int R) {
        return (((long) L) << 32) ^ (R & 0xffffffffL);
    }
}
