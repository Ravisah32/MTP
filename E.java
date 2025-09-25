import java.io.*;
import java.util.*;

public class E {

    // ----- Point -----
    static class Point {
        double[] coords = new double[3]; // rb, rg, bg
        int index;
    }

    // ----- Range Tree Node -----
    static class RangeTreeNode {
        int dim; // splitting dimension (0=rb, 1=rg, 2=bg)
        double split;
        Point pt;
        RangeTreeNode left, right;

        RangeTreeNode(int d, double s, Point p) {
            this.dim = d;
            this.split = s;
            this.pt = p;
        }
    }

    // Build range tree recursively
    static RangeTreeNode buildTree(List<Point> pts, int depth) {
        if (pts.isEmpty())
            return null;
        int dim = depth % 3;
        pts.sort(Comparator.comparingDouble(a -> a.coords[dim]));
        int mid = pts.size() / 2;
        RangeTreeNode node = new RangeTreeNode(dim, pts.get(mid).coords[dim], pts.get(mid));
        node.left = buildTree(new ArrayList<>(pts.subList(0, mid)), depth + 1);
        node.right = buildTree(new ArrayList<>(pts.subList(mid + 1, pts.size())), depth + 1);
        return node;
    }

    // Check if point inside query box
    static boolean inBox(Point p, double[][] box) {
        for (int i = 0; i < 3; i++) {
            if (p.coords[i] < box[i][0] || p.coords[i] > box[i][1])
                return false;
        }
        return true;
    }

    // Range query
    static void rangeQuery(RangeTreeNode root, double[][] box, List<Integer> result, int depth) {
        if (root == null)
            return;
        Point p = root.pt;
        if (inBox(p, box))
            result.add(p.index);

        int dim = root.dim;
        if (box[dim][0] <= root.split)
            rangeQuery(root.left, box, result, depth + 1);
        if (box[dim][1] >= root.split)
            rangeQuery(root.right, box, result, depth + 1);
    }

    // Trim whitespace
    static String trim(String s) {
        return s == null ? "" : s.trim();
    }

    // Compute Jaccard similarity between two ranges
    static double jaccard(int a1, int a2, int b1, int b2) {
        int interLeft = Math.max(a1, b1);
        int interRight = Math.min(a2, b2);
        int intersection = Math.max(0, interRight - interLeft + 1);

        int lenA = a2 - a1 + 1;
        int lenB = b2 - b1 + 1;
        int uni = lenA + lenB - intersection;

        if (uni == 0)
            return 0.0;
        return (double) intersection / uni;
    }

    public static void main(String[] args) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader("output_colors.csv"));
        String header = br.readLine(); // skip header

        List<Point> points = new ArrayList<>();
        String line;
        while ((line = br.readLine()) != null) {
            if (line.isEmpty())
                continue;
            String[] parts = line.split(",");
            if (parts.length < 5)
                continue;
            if (trim(parts[0]).isEmpty() || trim(parts[2]).isEmpty()
                    || trim(parts[3]).isEmpty() || trim(parts[4]).isEmpty())
                continue;

            Point p = new Point();
            p.index = Integer.parseInt(trim(parts[0]));
            p.coords[0] = Double.parseDouble(trim(parts[2])); // rb
            p.coords[1] = Double.parseDouble(trim(parts[3])); // rg
            p.coords[2] = Double.parseDouble(trim(parts[4])); // bg
            points.add(p);
        }
        br.close();

        // Build tree and measure time
        long buildStart = System.currentTimeMillis();
        RangeTreeNode root = buildTree(points, 0);
        long buildEnd = System.currentTimeMillis();
        System.out.println("Tree built in " + (buildEnd - buildStart) + " ms");

        // Map index -> Point
        Map<Integer, Point> indexToPoint = new HashMap<>();
        for (Point p : points) {
            indexToPoint.put(p.index, p);
        }

        Scanner sc = new Scanner(System.in);
        while (true) {
            System.out.println("Enter x and range [start end] (enter -1 -1 -1 to exit):");
            double x = sc.nextDouble();
            int start = sc.nextInt();
            int end = sc.nextInt();

            if (x == -1 && start == -1 && end == -1)
                break;

            if (start > end) {
                int tmp = start;
                start = end;
                end = tmp;
            }

            long queryStart = System.currentTimeMillis();

            double bestSim = -1.0;
            int bestQ = -1;
            int bestC = -1;

            // Iterate over all indices
            for (Point qPoint : points) {
                int qidx = qPoint.index;

                // build query center = coords of qidx-1 or 0,0,0
                double[] center = new double[] { 0.0, 0.0, 0.0 };
                if (qidx > 1) {
                    Point prev = indexToPoint.get(qidx - 1);
                    if (prev != null) {
                        center[0] = prev.coords[0];
                        center[1] = prev.coords[1];
                        center[2] = prev.coords[2];
                    }
                }

                double[][] box = new double[3][2];
                for (int i = 0; i < 3; i++) {
                    box[i][0] = center[i] - x;
                    box[i][1] = center[i] + x;
                }

                List<Integer> valid = new ArrayList<>();
                rangeQuery(root, box, valid, 0);

                // Exclude self, only > qidx
                valid.removeIf(v -> v == qidx || v <= qidx);
                Collections.sort(valid);

                for (int c : valid) {
                    int altStart = Math.min(qidx, c);
                    int altEnd = Math.max(qidx, c);
                    double sim = jaccard(start, end, altStart, altEnd);

                    if (sim > bestSim) {
                        bestSim = sim;
                        bestQ = qidx;
                        bestC = c;
                    }
                }
            }

            long queryEnd = System.currentTimeMillis();

            if (bestSim >= 0) {
                System.out.println("Best range = [" + bestQ + "," + bestC + "] (similarity = " + bestSim + ")");
            } else {
                System.out.println("No valid candidate found.");
            }

            System.out.println("Query executed in " + (queryEnd - queryStart) + " ms");
        }
        sc.close();
    }
}
