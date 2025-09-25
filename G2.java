import java.io.*;
import java.util.*;

public class G2 {

    public static void main(String[] args) throws IOException {
        if (args.length < 5) {
            System.out.println(
                    "Usage: java G1 <csvFile> <startMin> <startMax> <numTypes> <epsilon>");
            return;
        }

        String csvFile = args[0];
        double startMin = Double.parseDouble(args[1]);
        double startMax = Double.parseDouble(args[2]);
        int numTypes = Integer.parseInt(args[3]);
        int epsilon = Integer.parseInt(args[4]);

        // Try to read a dataset that already contains per-color cumulative counts
        // ("*Count" columns).
        // If not found, fall back to the legacy reader (salary,category) and build
        // prefix in-memory.
        Data data = readCSVCounts(csvFile);
        if (data == null) {
            // legacy fallback
            List<Person> people = readCSVLegacy(csvFile);
            // sort by salary
            people.sort(Comparator.comparingDouble(p -> p.salary));
            double[] salaries = new double[people.size()];
            for (int i = 0; i < people.size(); i++)
                salaries[i] = people.get(i).salary;

            // map categories to index
            Set<String> categorySet = new LinkedHashSet<>();
            for (Person p : people)
                categorySet.add(p.category);
            List<String> categories = new ArrayList<>(categorySet);
            Map<String, Integer> catIndex = new HashMap<>();
            for (int i = 0; i < categories.size(); i++)
                catIndex.put(categories.get(i), i);

            int C = categories.size();
            int N = people.size();
            int[][] prefix = new int[C][N + 1];
            for (int i = 0; i < N; i++) {
                int c = catIndex.get(people.get(i).category);
                for (int j = 0; j < C; j++) {
                    prefix[j][i + 1] = prefix[j][i] + (j == c ? 1 : 0);
                }
            }
            data = new Data(salaries, categories, prefix);
        }

        Range startRange = new Range(startMin, startMax);

        System.out.println("Counts in original range:");
        printCategoryCounts(data.salaries, data.prefix, data.categories, startRange);

        long queryStart = System.currentTimeMillis();
        Range fairRange = bfsFindFirstFairRange(data.salaries, data.prefix, data.categories, startRange, numTypes,
                epsilon);
        long queryEnd = System.currentTimeMillis();
        System.out.println("Query executed in " + (queryEnd - queryStart) + " ms");

        if (fairRange != null) {
            System.out.println("Fair range found: [" + fairRange.min + ", " + fairRange.max + "]");
            System.out.println("Counts in final fair range:");
            printCategoryCounts(data.salaries, data.prefix, data.categories, fairRange);
        } else {
            System.out.println("No fair range found.");
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

        Range(double a, double b) {
            this.min = a;
            this.max = b;
        }
    }

    static class Node {
        int L, R; // inclusive indices in salaries[]
        double sim; // priority score

        Node(int L, int R, double sim) {
            this.L = L;
            this.R = R;
            this.sim = sim;
        }
    }

    static class Data {
        final double[] salaries; // monotone non-decreasing scale (e.g., index or salary)
        final List<String> categories; // color names (derived from *Count headers)
        final int[][] prefix; // prefix[c][i] cumulative; we use [C][N+1]

        Data(double[] s, List<String> cats, int[][] p) {
            this.salaries = s;
            this.categories = cats;
            this.prefix = p;
        }
    }

    // ======== CSV reading (counts-first) ========

    // Reads CSV files that already include per-color cumulative counts in columns
    // named like "RedCount", "BlueCount", ...
    // Expected columns (order can vary):
    // - an index-like numeric column ("index" is common; if missing we use row
    // number 1..N)
    // - optional "color" column (ignored here)
    // - one or more *Count columns (cumulative)
    static Data readCSVCounts(String filename) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
            String header = br.readLine();
            if (header == null)
                return null;
            String[] cols = splitCSV(header);
            List<Integer> countColIdx = new ArrayList<>();
            List<String> colorNames = new ArrayList<>();
            Integer indexCol = null; // prefer a column literally named "index"; else infer later

            for (int i = 0; i < cols.length; i++) {
                String name = cols[i].trim();
                if (name.equalsIgnoreCase("index"))
                    indexCol = i;
                if (name.endsWith("Count")) {
                    countColIdx.add(i);
                    colorNames.add(name.substring(0, name.length() - "Count".length()));
                }
            }
            if (countColIdx.isEmpty())
                return null; // not a counts file; let legacy reader handle it

            // Read all rows
            List<double[]> rowsCounts = new ArrayList<>();
            List<Double> indexValues = new ArrayList<>();
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty())
                    continue;
                String[] parts = splitCSV(line);
                // Index value
                double idxVal;
                if (indexCol != null && indexCol < parts.length) {
                    idxVal = parseDoubleSafe(parts[indexCol]);
                } else {
                    // fallback: 1-based row number if no explicit index
                    idxVal = indexValues.size() + 1;
                }
                indexValues.add(idxVal);

                double[] counts = new double[countColIdx.size()];
                for (int k = 0; k < countColIdx.size(); k++) {
                    int cIdx = countColIdx.get(k);
                    counts[k] = (cIdx < parts.length) ? parseDoubleSafe(parts[cIdx]) : 0.0;
                }
                rowsCounts.add(counts);
            }

            int N = indexValues.size();
            int C = countColIdx.size();
            if (N == 0 || C == 0)
                return null;

            // salaries[] is the sortable/key axis. Here we use the index column values.
            double[] salaries = new double[N];
            for (int i = 0; i < N; i++)
                salaries[i] = indexValues.get(i);

            // Build prefix with size [C][N+1] such that prefix[c][i+1] == cumulative count
            // at row i
            int[][] prefix = new int[C][N + 1];
            for (int i = 0; i < N; i++) {
                double[] row = rowsCounts.get(i);
                for (int c = 0; c < C; c++) {
                    int cum = (int) Math.round(row[c]);
                    prefix[c][i + 1] = cum; // dataset already cumulative
                }
            }

            // Ensure color names are stable (empty names -> generic)
            for (int i = 0; i < colorNames.size(); i++) {
                String nm = colorNames.get(i);
                if (nm == null || nm.isEmpty())
                    colorNames.set(i, "Color" + (i + 1));
            }

            return new Data(salaries, colorNames, prefix);
        } catch (FileNotFoundException e) {
            throw e;
        }
    }

    // Legacy reader: salary,category (no prefix columns in file) -> we compute
    // prefix here
    static List<Person> readCSVLegacy(String filename) throws IOException {
        List<Person> people = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
            String line = br.readLine(); // header
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty())
                    continue;
                String[] parts = splitCSV(line);
                if (parts.length < 2)
                    continue;
                double sal = parseDoubleSafe(parts[0]);
                String cat = parts[1].trim();
                people.add(new Person(sal, cat));
            }
        }
        return people;
    }

    static double parseDoubleSafe(String s) {
        try {
            return Double.parseDouble(s.trim());
        } catch (Exception e) {
            return 0.0;
        }
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

    // ======== BFS / Informed BFS search ========

    static Range bfsFindFirstFairRange(
            double[] salaries,
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

        final int origSize = (R0 - L0 + 1);

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

        boolean[] seen = new boolean[n * 2 + 5]; // not used; we use a HashSet on packed keys
        HashSet<Long> visited = new HashSet<>();

        double sim0 = jaccardByIndex(L0, R0, L0, R0);
        pq.add(new Node(L0, R0, sim0));
        visited.add(packKey(L0, R0));

        while (!pq.isEmpty()) {
            Node cur = pq.poll();
            int L = cur.L, R = cur.R;
            if (isFair(salaries, prefix, L, R, categories.size(), numTypes, epsilon)) {
                return new Range(salaries[L], salaries[R]);
            }

            for (int[] nb : neighborsByIndex(salaries, L, R)) {
                int nL = nb[0], nR = nb[1];
                if (nL < 0 || nR >= n || nL > nR)
                    continue;
                long key = packKey(nL, nR);
                if (visited.contains(key))
                    continue;
                visited.add(key);
                double sim = jaccardByIndex(L0, R0, nL, nR);
                pq.add(new Node(nL, nR, sim));
            }
        }
        return null;
    }

    // Fairness: at least 'numTypes' present in [L..R], and max-min <= epsilon
    static boolean isFair(double[] salaries, int[][] prefix, int L, int R, int C, int numTypes, int epsilon) {
        int present = 0;
        int minCnt = Integer.MAX_VALUE, maxCnt = Integer.MIN_VALUE;
        for (int c = 0; c < C; c++) {
            int cnt = prefix[c][R + 1] - prefix[c][L];
            if (cnt > 0)
                present++;
            if (cnt < minCnt)
                minCnt = cnt;
            if (cnt > maxCnt)
                maxCnt = cnt;
        }
        if (present < numTypes)
            return false;
        if (maxCnt == 0)
            return false; // empty interval safeguard
        return (maxCnt - minCnt) <= epsilon;
    }

    // ======== Neighbor generation on index bounds ========
    // Returns a list of [L', R'] neighbors by minimally moving a boundary to the
    // next distinct key
    static List<int[]> neighborsByIndex(double[] s, int L, int R) {
        int n = s.length;
        List<int[]> v = new ArrayList<>(4);

        int nL = nextRightDistinctLeft(s, L, R);
        if (nL != -1)
            v.add(new int[] { nL, R });

        int pL = prevLeftDistinctLeft(s, L);
        if (pL != -1)
            v.add(new int[] { pL, R });

        int nR = nextRightDistinctRight(s, L, R);
        if (nR != -1)
            v.add(new int[] { L, nR });

        int pR = prevLeftDistinctRight(s, R);
        if (pR != -1)
            v.add(new int[] { L, pR });

        return v;
    }

    static int nextRightDistinctLeft(double[] s, int L, int R) {
        int n = s.length;
        if (L >= R)
            return -1;
        double val = s[L];
        int i = L + 1;
        while (i <= R && s[i] == val)
            i++;
        return (i <= R) ? i : -1;
    }

    static int prevLeftDistinctLeft(double[] s, int L) {
        if (L <= 0)
            return -1;
        double val = s[L];
        int i = L - 1;
        while (i >= 0 && s[i] == val)
            i--;
        return (i >= 0) ? i : -1;
    }

    static int nextRightDistinctRight(double[] s, int L, int R) {
        int n = s.length;
        if (L >= R)
            return -1;
        double val = s[R];
        int i = R - 1;
        while (i >= L && s[i] == val)
            i--;
        return (i >= L) ? i : -1;
    }

    static int prevLeftDistinctRight(double[] s, int R) {
        int n = s.length;
        if (R >= n - 1)
            return -1;
        double val = s[R];
        int i = R + 1;
        while (i < n && s[i] == val)
            i++;
        return (i < n) ? i : -1;
    }

    // ======== Jaccard similarity on index intervals ========

    static double jaccardByIndex(int aL, int aR, int bL, int bR) {
        int inter = intervalIntersectionSize(aL, aR, bL, bR);
        int ua = (aR - aL + 1) + (bR - bL + 1) - inter;
        return (ua == 0) ? 0.0 : ((double) inter) / ua;
    }

    static int intervalIntersectionSize(int aL, int aR, int bL, int bR) {
        int L = Math.max(aL, bL);
        int R = Math.min(aR, bR);
        return (L <= R) ? (R - L + 1) : 0;
    }

    static int lowerBound(double[] a, double x) {
        int n = a.length;
        int lo = 0, hi = n;
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            if (a[mid] < x)
                lo = mid + 1;
            else
                hi = mid;
        }
        return lo;
    }

    static int upperBound(double[] a, double x) {
        int n = a.length;
        int lo = 0, hi = n;
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            if (a[mid] <= x)
                lo = mid + 1;
            else
                hi = mid;
        }
        return lo;
    }

    static long packKey(int L, int R) {
        return (((long) L) << 32) ^ (R & 0xffffffffL);
    }
}
