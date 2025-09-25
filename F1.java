import java.io.*;
import java.util.*;

/**
 * E.java — brute-force search version that supports ANY number of colors.
 *
 * Input (same as G1/F1):
 * java E <csvFile> <startMin> <startMax> <numTypes> <epsilon>
 *
 * CSV expectations:
 * - Dataset already contains cumulative prefix columns per color named like
 * "RedCount", "BlueCount", ...
 * - Optionally an "index" column; if missing, row index 1..N is used as the key
 * axis.
 *
 * Fairness check:
 * - For any interval [L,R], compute per-color counts using the cumulative
 * columns:
 * cnt_c(L,R) = pref_c[R] - (L > 0 ? pref_c[L-1] : 0)
 * - Let maxCnt = max_c cnt_c(L,R), minCnt = min_c cnt_c(L,R)
 * - Interval is fair iff (maxCnt - minCnt) <= epsilon AND at least numTypes
 * colors have cnt > 0.
 *
 * Selection (ties resolved deterministically):
 * - Maximize Jaccard overlap with the start interval [L0,R0]
 * - Then minimize width (R-L)
 * - Then lexicographically by (L,R)
 *
 * Timing: program-wide in milliseconds (from start of main until the fair range
 * is computed).
 */
public class F1 {

    public static void main(String[] args) throws Exception {
        // --- Program-wide timer starts immediately ---
        long tProgramStart = System.nanoTime();

        if (args.length < 5) {
            System.out.println("Usage: java E <csvFile> <startMin> <startMax> <numTypes> <epsilon>");
            return;
        }

        final String csvFile = args[0];
        final double startMin = Double.parseDouble(args[1]);
        final double startMax = Double.parseDouble(args[2]);
        final int numTypes = Integer.parseInt(args[3]);
        final int epsilon = Integer.parseInt(args[4]);

        Data data = readCSVCounts(csvFile);
        if (data == null) {
            throw new IllegalArgumentException(
                    "CSV must contain cumulative *Count columns (e.g., RedCount, BlueCount, ...)");
        }

        Range startRange = new Range(startMin, startMax);
        System.out.println("Counts in original range:");
        printCategoryCounts(data.salaries, data.prefix, data.categories, startRange);

        // --- Run brute-force search ---
        Range fair = bruteForceBestFairRange(data, startRange, numTypes, epsilon);

        // --- Stop timer as soon as the result is computed (includes CSV read + search)
        // ---
        long tProgramEnd = System.nanoTime();
        double totalMs = (tProgramEnd - tProgramStart) / 1_000_000.0;
        System.out.printf("Total time (program): %.3f ms%n", totalMs);

        if (fair != null) {
            System.out.println("Fair range found: [" + fair.min + ", " + fair.max + "]");
            System.out.println("Counts in final fair range:");
            printCategoryCounts(data.salaries, data.prefix, data.categories, fair);
        } else {
            System.out.println("No fair range found.");
        }
    }

    // ================= data holders =================

    static final class Range {
        final double min, max;

        Range(double a, double b) {
            this.min = a;
            this.max = b;
        }
    }

    static final class Data {
        final double[] salaries; // index/key column (monotone, used for bounds)
        final List<String> categories; // names inferred from *Count columns
        final int[][] prefix; // prefix[c][i+1] == cumulative at row i

        Data(double[] s, List<String> cats, int[][] p) {
            this.salaries = s;
            this.categories = cats;
            this.prefix = p;
        }
    }

    // ================= brute force core =================

    static Range bruteForceBestFairRange(Data data, Range startRange, int numTypes, int epsilon) {
        double[] s = data.salaries;
        int[][] pref = data.prefix;
        int C = data.categories.size();
        int n = s.length;
        if (n == 0)
            return null;

        int L0 = lowerBound(s, startRange.min);
        int R0 = upperBound(s, startRange.max) - 1;
        if (L0 < 0)
            L0 = 0;
        if (R0 >= n)
            R0 = n - 1;
        if (L0 > R0)
            return null;

        double bestScore = -1.0;
        int bestL = -1, bestR = -1;

        for (int L = 0; L < n; L++) {
            for (int R = L; R < n; R++) {
                if (!isFair(pref, L, R, C, numTypes, epsilon))
                    continue;
                double score = jaccardByIndex(L0, R0, L, R);
                if (score > bestScore) {
                    bestScore = score;
                    bestL = L;
                    bestR = R;
                } else if (score == bestScore && bestScore >= 0) {
                    int wBest = bestR - bestL, wCur = R - L;
                    if (wCur < wBest || (wCur == wBest && (L < bestL || (L == bestL && R < bestR)))) {
                        bestL = L;
                        bestR = R;
                    }
                }
            }
        }

        if (bestL == -1)
            return null;
        return new Range(s[bestL], s[bestR]);
    }

    static boolean isFair(int[][] prefix, int L, int R, int C, int numTypes, int epsilon) {
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
        if (maxCnt == Integer.MIN_VALUE)
            return false; // safety
        return (maxCnt - minCnt) <= epsilon;
    }

    // ================= CSV reading =================

    static Data readCSVCounts(String filename) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
            String header = br.readLine();
            if (header == null)
                return null;
            String[] cols = splitCSV(header);

            List<Integer> countColIdx = new ArrayList<>();
            List<String> colorNames = new ArrayList<>();
            Integer indexCol = null;

            for (int i = 0; i < cols.length; i++) {
                String name = cols[i].trim();
                if (name.equalsIgnoreCase("index"))
                    indexCol = i;
                if (name.endsWith("Count")) {
                    countColIdx.add(i);
                    String base = name.substring(0, name.length() - "Count".length());
                    colorNames.add(base.isEmpty() ? ("Color" + (countColIdx.size())) : base);
                }
            }
            if (countColIdx.isEmpty())
                return null;

            List<Double> key = new ArrayList<>();
            List<int[]> rows = new ArrayList<>();
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty())
                    continue;
                String[] parts = splitCSV(line);

                double idxVal = (indexCol != null && indexCol < parts.length) ? parseDoubleSafe(parts[indexCol])
                        : (key.size() + 1);
                key.add(idxVal);

                int[] cum = new int[countColIdx.size()];
                for (int k = 0; k < countColIdx.size(); k++) {
                    int cIdx = countColIdx.get(k);
                    int v = 0;
                    if (cIdx < parts.length && parts[cIdx] != null && !parts[cIdx].isEmpty()) {
                        double d = parseDoubleSafe(parts[cIdx]);
                        v = (int) Math.round(d);
                    }
                    cum[k] = v;
                }
                rows.add(cum);
            }

            int n = key.size();
            int C = countColIdx.size();
            if (n == 0)
                return new Data(new double[0], colorNames, new int[C][1]);

            double[] salaries = new double[n];
            for (int i = 0; i < n; i++)
                salaries[i] = key.get(i);

            int[][] prefix = new int[C][n + 1];
            for (int i = 0; i < n; i++) {
                int[] cum = rows.get(i);
                for (int c = 0; c < C; c++) {
                    prefix[c][i + 1] = cum[c]; // file already cumulative
                }
            }

            return new Data(salaries, colorNames, prefix);
        }
    }

    static double parseDoubleSafe(String s) {
        try {
            return Double.parseDouble(s.trim());
        } catch (Exception e) {
            return 0.0;
        }
    }

    // simple CSV splitter (handles quoted commas)
    static String[] splitCSV(String line) {
        List<String> out = new ArrayList<>();
        StringBuilder sb = new StringBuilder();
        boolean inQ = false;
        for (int i = 0; i < line.length(); i++) {
            char ch = line.charAt(i);
            if (ch == '"') {
                inQ = !inQ;
            } else if (ch == ',' && !inQ) {
                out.add(sb.toString());
                sb.setLength(0);
            } else {
                sb.append(ch);
            }
        }
        out.add(sb.toString());
        return out.toArray(new String[0]);
    }

    // ================= utilities =================

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

    static int lowerBound(double[] a, double x) {
        int lo = 0, hi = a.length;
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
        int lo = 0, hi = a.length;
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            if (a[mid] <= x)
                lo = mid + 1;
            else
                hi = mid;
        }
        return lo;
    }

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
}
