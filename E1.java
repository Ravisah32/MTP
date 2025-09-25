import java.io.*;
import java.util.*;

/**
 * E.java — Interactive range search using an axis-aligned "range tree" over
 * m-dimensional prefix vectors built from pairwise cumulative columns (e.g.,
 * "Red-Blue").
 *
 * Usage (takes input filename from user, as requested):
 * java E <input.csv> <number_of_colors>
 *
 * Program flow:
 * 1) Read CSV; detect pairwise columns named like "A-B". Infer color set and
 * validate its size
 * equals <number_of_colors>.
 * 2) Build prefix vectors V_t in R^m for t = 0..n where V_0 = 0 and V_t holds
 * cumulative values
 * for each pairwise column at row t.
 * 3) Build an axis-aligned range tree over the points { (V_t, t) }.
 * 4) Print: "Tree built in <ms>"
 * 5) Enter interactive loop with the *exact* prompt you showed:
 * Enter x and range [start end] (enter -1 -1 -1 to exit):
 * Here x = epsilon (fairness threshold), start/end are *value* endpoints on the
 * key axis (index column).
 * 6) For each query, return the fair interval [L,R] with maximum Jaccard vs the
 * mapped [L0,R0].
 * Print exactly:
 * Best range = [<leftVal>,<rightVal>] (similarity = <val>)
 * Query executed in <ms> ms
 */
public class E1 {

    public static void main(String[] args) throws Exception {
        if (args.length < 2) {
            System.out.println("Usage: java E <input.csv> <number_of_colors>");
            return;
        }

        final String csvFile = args[0];
        final int numberOfColors = Integer.parseInt(args[1]);

        long t0 = System.nanoTime();
        Data data = readCSV(csvFile);
        if (data.m == 0) {
            System.out.println("No pairwise columns like 'Red-Blue' found in header.");
            return;
        }
        if (data.uniqueColors.size() != numberOfColors) {
            System.out.println("number_of_colors mismatch: file has " + data.uniqueColors.size() +
                    " unique colors but arg was " + numberOfColors);
            return;
        }

        RangeTree tree = new RangeTree(data.points, data.m); // points include V_0..V_n
        long t1 = System.nanoTime();
        System.out.println("Tree built in " + ((t1 - t0) / 1_000_000) + " ms");

        // Interactive query loop
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        while (true) {
            System.out.print("Enter x and range [start end] (enter -1 -1 -1 to exit):\n");
            String line = br.readLine();
            if (line == null)
                break;
            line = line.trim();
            if (line.isEmpty())
                continue;
            String[] parts = line.split("\\s+");
            if (parts.length < 3)
                continue;
            int epsilon, iStart = 0, iEnd = 0;
            try {
                epsilon = Integer.parseInt(parts[0]);
                iStart = (parts.length >= 2) ? Integer.parseInt(parts[1]) : 0;
                iEnd = (parts.length >= 3) ? Integer.parseInt(parts[2]) : 0;
            } catch (Exception ex) {
                continue;
            }
            if (epsilon == -1 && iStart == -1 && iEnd == -1)
                break;

            double startVal = iStart;
            double endVal = iEnd;

            long q0 = System.nanoTime();
            Result res = queryBest(data, tree, epsilon, startVal, endVal);
            long q1 = System.nanoTime();

            if (res == null) {
                System.out.println("Best range = [NA,NA] (similarity = 0.0)");
                System.out.println("Query executed in " + ((q1 - q0) / 1_000_000) + " ms");
                continue;
            }
            System.out.println(
                    "Best range = [" + res.leftVal + "," + res.rightVal + "] (similarity = " + res.similarity + ")");
            System.out.println("Query executed in " + ((q1 - q0) / 1_000_000) + " ms");
        }
    }

    // ===================== Query logic =====================

    static final class Result {
        final double leftVal, rightVal;
        final double similarity;

        Result(double l, double r, double s) {
            leftVal = l;
            rightVal = r;
            similarity = s;
        }
    }

    static Result queryBest(Data data, RangeTree tree, int epsilon, double startVal, double endVal) {
        int n = data.n;
        if (n == 0)
            return null;

        // Map value endpoints to index interval [L0,R0] (inclusive)
        int L0 = lowerBound(data.keys, Math.min(startVal, endVal));
        int R0 = upperBound(data.keys, Math.max(startVal, endVal)) - 1;
        if (L0 < 0)
            L0 = 0;
        if (R0 >= n)
            R0 = n - 1;
        if (L0 > R0)
            return null;

        double bestJ = -1.0;
        int bestL = -1, bestR = -1;

        // Iterate left boundary via prefix index p = L-1 in 0..n-1
        for (int p = 0; p <= n - 1; p++) {
            double[] center = data.points[p];
            double[] lo = new double[data.m];
            double[] hi = new double[data.m];
            for (int d = 0; d < data.m; d++) {
                lo[d] = center[d] - epsilon;
                hi[d] = center[d] + epsilon;
            }

            // Range query returns indices t in [0..n], where t = R (prefix index)
            List<Integer> hits = new ArrayList<>();
            tree.rangeQuery(lo, hi, hits);

            // For each candidate R >= p+1 (strictly to the right), evaluate Jaccard on
            // *index* axis
            for (int id : hits) {
                int R = id;
                if (R <= p)
                    continue; // ensure L <= R
                int L = p + 1;
                double j = jaccard(L, R, L0 + 1, R0 + 1); // use 1-based lengths for clarity
                if (j > bestJ) {
                    bestJ = j;
                    bestL = L;
                    bestR = R;
                } else if (j == bestJ && j >= 0) {
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
        // Map back to value endpoints (data.keys are 0-based rows; bestL..bestR are
        // 1-based positions)
        double leftVal = data.keys[bestL - 1];
        double rightVal = data.keys[bestR - 1];
        return new Result(leftVal, rightVal, bestJ);
    }

    // ===================== Data loading =====================

    static final class Data {
        final int n; // number of rows
        final double[] keys; // key axis values (length n)
        final int m; // dimensions = number of pair columns
        final List<String> pairNames;// e.g., "Red-Blue"
        final Set<String> uniqueColors;
        final int[][] prefix; // [m][n+1]
        final double[][] points; // V_t in R^m for t=0..n (points[0] = 0^m)

        Data(int n, double[] keys, int m, List<String> pairNames, Set<String> uniqueColors,
                int[][] prefix, double[][] points) {
            this.n = n;
            this.keys = keys;
            this.m = m;
            this.pairNames = pairNames;
            this.uniqueColors = uniqueColors;
            this.prefix = prefix;
            this.points = points;
        }
    }

    static Data readCSV(String filename) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
            String header = br.readLine();
            if (header == null)
                throw new IllegalArgumentException("Empty CSV");
            String[] cols = splitCSV(header);

            Integer indexCol = null;
            List<Integer> pairCols = new ArrayList<>();
            List<String> pairNames = new ArrayList<>();
            Set<String> colors = new LinkedHashSet<>();

            for (int i = 0; i < cols.length; i++) {
                String name = cols[i].trim();
                if (name.equalsIgnoreCase("index"))
                    indexCol = i;
                else if (isPairName(name)) {
                    pairCols.add(i);
                    pairNames.add(name);
                    String[] ab = name.split("-");
                    if (ab.length == 2) {
                        colors.add(ab[0]);
                        colors.add(ab[1]);
                    }
                }
            }
            if (pairCols.isEmpty())
                return new Data(0, new double[0], 0, pairNames, colors, new int[0][0], new double[0][0]);

            List<double[]> rows = new ArrayList<>();
            List<Double> key = new ArrayList<>();

            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty())
                    continue;
                String[] parts = splitCSV(line);
                double idxVal;
                if (indexCol != null && indexCol < parts.length)
                    idxVal = parseDoubleSafe(parts[indexCol]);
                else {
                    Double maybe = tryParse(parts[0]);
                    idxVal = (maybe != null) ? maybe : (key.size() + 1);
                }
                key.add(idxVal);

                double[] vals = new double[pairCols.size()];
                for (int k = 0; k < pairCols.size(); k++) {
                    int cIdx = pairCols.get(k);
                    vals[k] = (cIdx < parts.length && !parts[cIdx].isEmpty()) ? parseDoubleSafe(parts[cIdx]) : 0.0;
                }
                rows.add(vals);
            }

            int n = rows.size();
            int m = pairCols.size();
            double[] keys = new double[n];
            for (int i = 0; i < n; i++)
                keys[i] = key.get(i);

            int[][] prefix = new int[m][n + 1];
            for (int t = 0; t < n; t++) {
                double[] row = rows.get(t);
                for (int d = 0; d < m; d++)
                    prefix[d][t + 1] = (int) Math.round(row[d]);
            }

            double[][] points = new double[n + 1][m]; // V_0..V_n
            for (int t = 1; t <= n; t++)
                for (int d = 0; d < m; d++)
                    points[t][d] = prefix[d][t];
            // points[0] stays 0s

            return new Data(n, keys, m, pairNames, colors, prefix, points);
        }
    }

    static boolean isPairName(String name) {
        int dash = name.indexOf('-');
        return dash > 0 && dash < name.length() - 1 && name.indexOf('-', dash + 1) == -1;
    }

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

    static Double tryParse(String s) {
        try {
            return Double.parseDouble(s.trim());
        } catch (Exception e) {
            return null;
        }
    }

    static double parseDoubleSafe(String s) {
        try {
            return Double.parseDouble(s.trim());
        } catch (Exception e) {
            return 0.0;
        }
    }

    // ===================== "Range tree" (axis-aligned, recursive split)
    // =====================
    // This structure supports orthogonal box queries in m-D. It recursively
    // partitions points by
    // median along the current dimension (cycling dims), and stores bounding boxes
    // for pruning.
    //
    // Although implementation resembles a kd-style recursive partition, we keep
    // strict axis-aligned
    // box querying semantics required for a range tree use case.

    static final class RangeTree {
        final int dims;
        final Node root;
        final double[][] P; // reference to points V_t

        RangeTree(double[][] points, int dims) {
            this.dims = dims;
            this.P = points;
            int n = points.length;
            int[] all = new int[n];
            for (int i = 0; i < n; i++)
                all[i] = i; // indices 0..n
            this.root = build(all, 0, n, 0);
        }

        final class Node {
            // bounding box
            double[] lo, hi;
            // split info
            int dim;
            double split;
            Node left, right;
            int[] idxs; // leaf indices
            boolean isLeaf;
        }

        Node build(int[] idx, int loi, int hii, int depth) {
            int n = hii - loi;
            Node node = new Node();
            node.lo = new double[dims];
            node.hi = new double[dims];
            Arrays.fill(node.lo, Double.POSITIVE_INFINITY);
            Arrays.fill(node.hi, Double.NEGATIVE_INFINITY);
            for (int k = loi; k < hii; k++) {
                double[] pt = P[idx[k]];
                for (int d = 0; d < dims; d++) {
                    if (pt[d] < node.lo[d])
                        node.lo[d] = pt[d];
                    if (pt[d] > node.hi[d])
                        node.hi[d] = pt[d];
                }
            }
            if (n <= 32) { // leaf threshold
                node.isLeaf = true;
                node.idxs = Arrays.copyOfRange(idx, loi, hii);
                return node;
            }
            int dim = depth % dims;
            int mid = loi + n / 2;
            nthElement(idx, loi, mid, hii, dim);
            node.dim = dim;
            node.split = P[idx[mid]][dim];
            node.left = build(idx, loi, mid, depth + 1);
            node.right = build(idx, mid, hii, depth + 1);
            return node;
        }

        void nthElement(int[] a, int lo, int mid, int hi, int dim) {
            int l = lo, r = hi - 1;
            while (true) {
                int i = l, j = r;
                double pivot = P[a[(l + r) >>> 1]][dim];
                while (i <= j) {
                    while (P[a[i]][dim] < pivot)
                        i++;
                    while (P[a[j]][dim] > pivot)
                        j--;
                    if (i <= j) {
                        int t = a[i];
                        a[i] = a[j];
                        a[j] = t;
                        i++;
                        j--;
                    }
                }
                if (j < mid)
                    l = i;
                else
                    r = j;
                if (l >= mid && r <= mid)
                    return;
            }
        }

        void rangeQuery(double[] lo, double[] hi, List<Integer> out) {
            rangeQuery(root, lo, hi, out);
        }

        void rangeQuery(Node node, double[] lo, double[] hi, List<Integer> out) {
            if (node == null)
                return;
            if (!overlaps(node, lo, hi))
                return; // prune by bounding box
            if (node.isLeaf) {
                for (int id : node.idxs)
                    if (contains(P[id], lo, hi))
                        out.add(id);
                return;
            }
            // Recurse both sides if needed; bounding boxes already prune a lot
            rangeQuery(node.left, lo, hi, out);
            rangeQuery(node.right, lo, hi, out);
        }

        boolean overlaps(Node node, double[] qlo, double[] qhi) {
            for (int d = 0; d < dims; d++)
                if (node.hi[d] < qlo[d] || node.lo[d] > qhi[d])
                    return false;
            return true;
        }

        boolean contains(double[] pt, double[] qlo, double[] qhi) {
            for (int d = 0; d < dims; d++)
                if (pt[d] < qlo[d] || pt[d] > qhi[d])
                    return false;
            return true;
        }
    }

    // ===================== Math & utils =====================

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

    static double jaccard(int aL, int aR, int bL, int bR) {
        int inter = intersectionSize(aL, aR, bL, bR);
        int uni = (aR - aL + 1) + (bR - bL + 1) - inter;
        return (uni == 0) ? 0.0 : ((double) inter) / uni;
    }

    static int intersectionSize(int aL, int aR, int bL, int bR) {
        int L = Math.max(aL, bL), R = Math.min(aR, bR);
        return (L <= R) ? (R - L + 1) : 0;
    }
}
