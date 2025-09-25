import java.io.*;
import java.util.*;

public class F {

    static class Row {
        int index;
        String color;
        int red, blue, green;
    }

    static String trim(String s) {
        return s == null ? "" : s.trim();
    }

    static int rangeIntersection(int a, int b, int c, int d) {
        return Math.max(0, Math.min(b, d) - Math.max(a, c) + 1);
    }

    static int rangeUnion(int a, int b, int c, int d) {
        int inter = rangeIntersection(a, b, c, d);
        return (b - a + 1) + (d - c + 1) - inter;
    }

    static class Range {
        int startIdx, endIdx;
        int red, blue, green;

        Range(int s, int e, int r, int b, int g) {
            startIdx = s;
            endIdx = e;
            red = r;
            blue = b;
            green = g;
        }
    }

    public static void main(String[] args) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader("brute_output.csv"));
        String line = br.readLine(); // skip header

        List<Row> rows = new ArrayList<>();
        while ((line = br.readLine()) != null) {
            if (line.isEmpty())
                continue;
            String[] parts = line.split(",");
            Row r = new Row();
            r.index = Integer.parseInt(trim(parts[0]));
            r.color = trim(parts[1]);
            r.red = Integer.parseInt(trim(parts[2]));
            r.blue = Integer.parseInt(trim(parts[3]));
            r.green = Integer.parseInt(trim(parts[4]));
            rows.add(r);
        }
        br.close();

        if (rows.isEmpty()) {
            System.err.println("No rows found in brute_output.csv");
            return;
        }

        Scanner sc = new Scanner(System.in);
        System.out.print("Enter x: ");
        int x = sc.nextInt();

        System.out.print("Enter input range (startIndex endIndex): ");
        int p = sc.nextInt(), q = sc.nextInt();
        if (p > q) {
            int tmp = p;
            p = q;
            q = tmp;
        }

        // ---- Measure query time ----
        long queryStart = System.currentTimeMillis();

        double bestSim = -1.0;
        List<Range> bestRanges = new ArrayList<>();

        int n = rows.size();
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                int redInRange = rows.get(j).red - rows.get(i).red +
                        (rows.get(i).color.equals("Red") ? 1 : 0);
                int blueInRange = rows.get(j).blue - rows.get(i).blue +
                        (rows.get(i).color.equals("Blue") ? 1 : 0);
                int greenInRange = rows.get(j).green - rows.get(i).green +
                        (rows.get(i).color.equals("Green") ? 1 : 0);

                if (Math.abs(redInRange - blueInRange) <= x &&
                        Math.abs(blueInRange - greenInRange) <= x &&
                        Math.abs(redInRange - greenInRange) <= x) {

                    int inter = rangeIntersection(p, q, rows.get(i).index, rows.get(j).index);
                    int uni = rangeUnion(p, q, rows.get(i).index, rows.get(j).index);
                    if (uni == 0)
                        continue;

                    double sim = (double) inter / uni;

                    if (sim > bestSim + 1e-9) {
                        bestSim = sim;
                        bestRanges.clear();
                        bestRanges.add(new Range(rows.get(i).index, rows.get(j).index,
                                redInRange, blueInRange, greenInRange));
                    } else if (Math.abs(sim - bestSim) < 1e-9) {
                        bestRanges.add(new Range(rows.get(i).index, rows.get(j).index,
                                redInRange, blueInRange, greenInRange));
                    }
                }
            }
        }

        long queryEnd = System.currentTimeMillis();

        if (bestRanges.isEmpty()) {
            System.out.println("No similar ranges found.");
        } else {
            System.out.println("Most similar ranges (Jaccard similarity = " + bestSim + "):");
            for (Range r : bestRanges) {
                System.out.println("[" + r.startIdx + "," + r.endIdx + "] " +
                        " Red=" + r.red + " Blue=" + r.blue + " Green=" + r.green);
            }
        }

        System.out.println("Query executed in " + (queryEnd - queryStart) + " ms");
    }
}
