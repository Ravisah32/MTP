import java.io.*;
import java.util.*;

public class G

{

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

            // build prefix counts
            int[][] prefix = new int[categories.size()][people.size() + 1];
            for (int i = 0; i < people.size(); i++) {
                int c = catIndex.get(people.get(i).category);
                for (int j = 0; j < categories.size(); j++) {
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

    static List<Person> readCSV(String filename) throws IOException {
        List<Person> people = new ArrayList<>();
        BufferedReader br = new BufferedReader(new FileReader(filename));
        br.readLine(); // skip header
        String line;
        while ((line = br.readLine()) != null) {
            line = line.trim();
            if (!line.isEmpty()) {
                String[] parts = line.split(",");
                if (parts.length >= 2) {
                    try {
                        double salary = Double.parseDouble(parts[0].trim());
                        String category = parts[1].trim();
                        people.add(new Person(salary, category));
                    } catch (NumberFormatException e) {
                        // skip invalid rows
                    }
                }
            }
        }
        br.close();
        return people;
    }

    static boolean isFair(double[] salaries, int[][] prefix, List<String> categories,
            Range r, int numTypes, int epsilon) {
        int L = Arrays.binarySearch(salaries, r.min);
        if (L < 0)
            L = -L - 1;
        int R = Arrays.binarySearch(salaries, r.max);
        if (R < 0)
            R = -R - 2;

        if (L > R)
            return false;

        int minCount = Integer.MAX_VALUE;
        int maxCount = Integer.MIN_VALUE;

        for (int c = 0; c < categories.size(); c++) {
            int cnt = prefix[c][R + 1] - prefix[c][L];
            minCount = Math.min(minCount, cnt);
            maxCount = Math.max(maxCount, cnt);
        }

        if (categories.size() != numTypes)
            return false;
        return maxCount - minCount <= epsilon;
    }

    static void printCategoryCounts(double[] salaries, int[][] prefix, List<String> categories, Range r) {
        int L = Arrays.binarySearch(salaries, r.min);
        if (L < 0)
            L = -L - 1;
        int R = Arrays.binarySearch(salaries, r.max);
        if (R < 0)
            R = -R - 2;

        if (L > R) {
            System.out.println("No salaries in range.");
            return;
        }

        for (int c = 0; c < categories.size(); c++) {
            int cnt = prefix[c][R + 1] - prefix[c][L];
            System.out.println(categories.get(c) + ": " + cnt);
        }
    }

    static Range bfsFindFirstFairRange(double[] salaries, int[][] prefix, List<String> categories,
            Range startRange, int numTypes, int epsilon) {
        LinkedList<Range> queue = new LinkedList<>();
        Set<String> visited = new HashSet<>();
        queue.add(startRange);
        visited.add(key(startRange));

        while (!queue.isEmpty()) {
            Range current = queue.poll();

            if (isFair(salaries, prefix, categories, current, numTypes, epsilon)) {
                return current; // ✅ stop immediately
            }

            int idxMin = findIndex(salaries, current.min, true);
            int idxMax = findIndex(salaries, current.max, false);

            List<Range> neighbors = new ArrayList<>();
            if (idxMin > 0) {
                neighbors.add(new Range(salaries[idxMin - 1], current.max));
            }
            if (idxMin < salaries.length - 1 && salaries[idxMin] < current.max) {
                neighbors.add(new Range(salaries[idxMin + 1], current.max));
            }
            if (idxMax < salaries.length - 1) {
                neighbors.add(new Range(current.min, salaries[idxMax + 1]));
            }
            if (idxMax > 0 && salaries[idxMax] > current.min) {
                neighbors.add(new Range(current.min, salaries[idxMax - 1]));
            }

            for (Range r : neighbors) {
                if (r.min < r.max) {
                    String k = key(r);
                    if (!visited.contains(k)) {
                        visited.add(k);
                        queue.add(r);
                    }
                }
            }
        }
        return null;
    }

    static String key(Range r) {
        return r.min + "-" + r.max;
    }

    static int findIndex(double[] salaries, double val, boolean left) {
        if (left) {
            for (int i = 0; i < salaries.length; i++) {
                if (salaries[i] >= val)
                    return i;
            }
        } else {
            for (int i = salaries.length - 1; i >= 0; i--) {
                if (salaries[i] <= val)
                    return i;
            }
        }
        return -1;
    }

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
}
