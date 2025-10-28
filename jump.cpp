#include <bits/stdc++.h>
using namespace std;

// ===================== Data structures =====================
class Point_1D
{
public:
    double X;       // 1-D coordinate (here: index from CSV)
    bool color;     // true = blue, false = red
    int cumulative; // prefix sum (+1 for blue, -1 for red)
    int forward_ptr, backward_ptr;

    Point_1D(double coordinate = 0.0, bool col = false)
        : X(coordinate), color(col), cumulative(0), forward_ptr(-1), backward_ptr(-1) {}

    bool operator<(const Point_1D &other) const { return X < other.X; }
};

class Range
{
public:
    int start, end; // indices into the sorted-with-sentinels vector
    Range(int st = 0, int e = 0) : start(st), end(e) {}
};

static inline int color_to_value(bool c) { return c ? 1 : -1; } // blue:+1, red:-1

// ===================== Preprocessing =====================
void preprocess(vector<Point_1D> *db)
{
    sort(db->begin(), db->end());

    // Sentinels
    db->insert(db->begin(), Point_1D(-numeric_limits<double>::infinity(), false));
    db->push_back(Point_1D(numeric_limits<double>::infinity(), false));

    // Backward pointers (right->left)
    map<int, vector<int>> ptr_map;
    int cumul = 0;
    for (int i = (int)db->size() - 1; i > 0; --i)
    {
        if ((int)db->size() - 1 == i)
        {
            db->at(i).cumulative = cumul;
            continue;
        }
        int cur = color_to_value(db->at(i).color);
        cumul += cur;
        db->at(i).cumulative = cumul;

        if (ptr_map.count(cumul - 2 * cur))
            ptr_map[cumul - 2 * cur].push_back(i + 1);
        else
            ptr_map.insert({cumul - 2 * cur, vector<int>{i + 1}});

        if (ptr_map.count(cumul))
        {
            for (int j : ptr_map[cumul])
                db->at(j).backward_ptr = i;
            ptr_map.erase(cumul);
        }
    }

    // Forward pointers (left->right)
    ptr_map.clear();
    cumul = 0;
    for (int i = 0; i < (int)db->size() - 1; ++i)
    {
        if (i == 0)
        {
            db->at(i).cumulative = cumul;
            continue;
        }
        int cur = color_to_value(db->at(i).color);
        cumul += cur;
        db->at(i).cumulative = cumul;

        if (ptr_map.count(cumul - 2 * cur))
            ptr_map[cumul - 2 * cur].push_back(i - 1);
        else
            ptr_map.insert({cumul - 2 * cur, vector<int>{i - 1}});

        if (ptr_map.count(cumul))
        {
            for (int j : ptr_map[cumul])
                db->at(j).forward_ptr = i;
            ptr_map.erase(cumul);
        }
    }
}

// ===================== Similarity =====================
static inline double Jaccard_similarity(Range r1, Range r2)
{
    int uni = max(r1.end, r2.end) - min(r1.start, r2.start) + 1;
    int inter = min(r1.end, r2.end) - max(r1.start, r2.start) + 1;
    if (uni <= 0)
        return 0.0;
    return double(max(0, inter)) / double(uni);
}

// ===================== Search helpers =====================
int binary_search_X(vector<Point_1D> *db, int start, int end, double x)
{
    if (end < start)
        return start;
    if (start == end)
        return start;
    int mid = (start + end) / 2;
    if (db->at(mid).X == x)
        return mid;
    if (db->at(mid).X < x)
        return binary_search_X(db, mid + 1, end, x);
    return binary_search_X(db, start, mid - 1, x);
}

// Your fast fair-range routine
Range get_fair_range(vector<Point_1D> *db, Range input_range, int eps,
                     double similarity(Range, Range))
{
    int disparity = db->at(input_range.end).cumulative - db->at(input_range.start - 1).cumulative;
    if (abs(disparity) <= eps)
        return input_range;

    double best_similarity = 0.0;
    Range fair_range = input_range;

    int disparity_difference = abs(disparity) - eps;
    int jumps_to_make = disparity_difference;

    stack<int> start_left, start_right, end_left, end_right;
    bool excess_color = disparity > 0; // true => too many blues

    start_left.push(input_range.start);
    start_right.push(input_range.start);
    end_left.push(input_range.end);
    end_right.push(input_range.end);

    // Advance start_right by required jumps
    while (jumps_to_make)
    {
        int top = start_right.top();
        if (db->at(top).color == excess_color)
            start_right.push(top + 1);
        else
            start_right.push(db->at(top - 1).forward_ptr + 1);
        --jumps_to_make;
    }

    // Explore right expansions of end
    jumps_to_make = disparity_difference;
    while (jumps_to_make)
    {
        int sr = start_right.top();
        int el = end_left.top();
        int er = end_right.top();

        Range sr_el(sr, el);
        double sim1 = similarity(input_range, sr_el);
        if (sim1 > best_similarity)
        {
            best_similarity = sim1;
            fair_range = sr_el;
        }

        if ((int)end_right.size() + (int)start_right.size() == disparity_difference + 2)
        {
            Range sr_er(sr, er);
            double sim2 = similarity(input_range, sr_er);
            if (sim2 > best_similarity)
            {
                best_similarity = sim2;
                fair_range = sr_er;
            }
        }

        if (er < (int)db->size() - 1)
        {
            if (db->at(er + 1).color != excess_color)
                end_right.push(er + 1);
            else if (db->at(er).forward_ptr != -1)
                end_right.push(db->at(er).forward_ptr);
        }
        if (db->at(el).color == excess_color)
            end_left.push(el - 1);
        else
            end_left.push(db->at(el + 1).backward_ptr - 1);

        start_right.pop();
        --jumps_to_make;
    }

    // Explore left expansions of start
    jumps_to_make = disparity_difference + 1;
    while (jumps_to_make)
    {
        int sl = start_left.top();
        int el = end_left.top();
        int er = end_right.top();

        Range sl_el(sl, el);
        double sim3 = similarity(input_range, sl_el);
        if (sim3 > best_similarity)
        {
            best_similarity = sim3;
            fair_range = sl_el;
        }

        if ((int)end_right.size() + (int)start_left.size() == disparity_difference + 2)
        {
            Range sl_er(sl, er);
            double sim4 = similarity(input_range, sl_er);
            if (sim4 > best_similarity)
            {
                best_similarity = sim4;
                fair_range = sl_er;
            }
            end_right.pop();
        }
        end_left.pop();

        if (sl > 1)
        {
            if (db->at(sl - 1).color != excess_color)
                start_left.push(sl - 1);
            else if (db->at(sl).backward_ptr != -1)
                start_left.push(db->at(sl).backward_ptr);
            else
                break;
        }
        else
            break;

        --jumps_to_make;
    }
    return fair_range;
}

// ===================== CSV/whitespace reader =====================
static vector<string> split_csv(const string &s, char delim = ',')
{
    vector<string> out;
    string cur;
    for (char c : s)
    {
        if (c == delim)
        {
            out.push_back(cur);
            cur.clear();
        }
        else
            cur.push_back(c);
    }
    out.push_back(cur);
    return out;
}
static bool parse_color_token(const string &tok, int &color01)
{
    string t = tok;
    for (auto &c : t)
        c = (char)tolower(c);
    if (t == "1" || t == "blue" || t == "true")
    {
        color01 = 1;
        return true;
    }
    if (t == "0" || t == "red" || t == "false")
    {
        color01 = 0;
        return true;
    }
    return false;
}

// Supports:
// 1) CSV header: "index,color,RedCount,BlueCount"  -> uses index as X
// 2) CSV: "color,value" (e.g., "blue,52345")
// 3) Space: "color value" or "blue 52345"
vector<Point_1D> *read_dataset(const string &filename)
{
    auto *database = new vector<Point_1D>();
    ifstream fh(filename);
    if (!fh.is_open())
    {
        cerr << "Error: failed to open file '" << filename << "'\n";
        return database;
    }
    string line;
    while (getline(fh, line))
    {
        if (line.empty())
            continue;

        // trim
        while (!line.empty() && isspace((unsigned char)line.back()))
            line.pop_back();
        size_t p = 0;
        while (p < line.size() && isspace((unsigned char)line[p]))
            ++p;
        if (p > 0)
            line = line.substr(p);

        int color01 = -1;
        double xval = 0.0;

        if (line.find(',') != string::npos)
        {
            auto toks = split_csv(line, ',');
            if (!toks.empty())
            {
                string t0 = toks[0];
                for (auto &c : t0)
                    c = (char)tolower(c);
                if (t0 == "index")
                    continue; // header line
            }

            if (toks.size() >= 4)
            {
                // index,color,RedCount,BlueCount
                try
                {
                    xval = stod(toks[0]);
                }
                catch (...)
                {
                    continue;
                }
                if (!parse_color_token(toks[1], color01))
                    continue;
            }
            else if (toks.size() == 2)
            {
                if (!parse_color_token(toks[0], color01))
                    continue;
                try
                {
                    xval = stod(toks[1]);
                }
                catch (...)
                {
                    continue;
                }
            }
            else
            {
                continue;
            }
        }
        else
        {
            istringstream iss(line);
            if (!(iss >> color01 >> xval))
            {
                iss.clear();
                iss.str(line);
                string cstr;
                double val;
                if ((iss >> cstr >> val) && parse_color_token(cstr, color01))
                    xval = val;
                else
                    continue;
            }
        }
        database->push_back(Point_1D(xval, color01 == 1)); // 1->blue, 0->red
    }
    return database;
}

// ===================== MAIN =====================
int main(int argc, char **argv)
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string fname;
    if (argc >= 2)
    {
        fname = argv[1];
    }
    else
    {
        cout << "Enter input filename: ";
        if (!getline(cin, fname) || fname.empty())
        {
            cerr << "No filename provided.\n";
            return 1;
        }
    }

    // Load data
    vector<Point_1D> *db = read_dataset(fname);
    if (!db || db->empty())
    {
        cerr << "No data loaded from '" << fname << "'.\n";
        delete db;
        return 1;
    }

    // Preprocess (NOT timed)
    preprocess(db);

    // Query input (coordinate range; for your CSV this means index range)
    double start_val, end_val;
    cout << "Enter query range (start end): ";
    if (!(cin >> start_val >> end_val))
    {
        cerr << "Invalid range input.\n";
        delete db;
        return 1;
    }
    if (start_val > end_val)
        swap(start_val, end_val);

    // Map to index range in the (sorted-with-sentinels) array
    int start_idx = binary_search_X(db, 1, (int)db->size() - 1, start_val);
    int end_idx = binary_search_X(db, 1, (int)db->size() - 1, end_val);
    start_idx = max(1, min(start_idx, (int)db->size() - 1));
    end_idx = max(1, min(end_idx, (int)db->size() - 1));
    if (start_idx > end_idx)
        swap(start_idx, end_idx);

    int eps = 0; // tolerance; change if you need slack

    // ---- TIME ONLY THE QUERY ----
    auto t0 = std::chrono::steady_clock::now();
    Range out = get_fair_range(db, Range(start_idx, end_idx), eps, Jaccard_similarity);
    auto t1 = std::chrono::steady_clock::now();
    double query_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // Report
    int in_disp = abs(db->at(end_idx).cumulative - db->at(start_idx - 1).cumulative);
    int out_disp = abs(db->at(out.end).cumulative - db->at(out.start - 1).cumulative);

    cout << fixed << setprecision(6);
    cout << "Query time (ms): " << query_ms << "\n";
    cout << "epsilon: " << eps << "\n";
    cout << "input disparity: " << in_disp << "\n";
    cout << "input indices: (" << start_idx << ", " << end_idx << ")\n";
    cout << "input range:   (" << db->at(start_idx).X << ", " << db->at(end_idx).X << ")\n";
    cout << "output indices:(" << out.start << ", " << out.end << ")\n";
    cout << "output range:  (" << db->at(out.start).X << ", " << db->at(out.end).X << ")\n";
    cout << "output disparity: " << out_disp << "\n";
    cout << "points in output: " << (out.end - out.start + 1) << "\n";

    delete db;
    return 0;
}
