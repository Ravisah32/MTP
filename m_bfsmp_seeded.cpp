#include <iostream>
#include <vector>
#include <limits>
#include <stack>
#include <algorithm>
#include <unordered_set>
#include <fstream>
#include <queue>
#include <bits/stdc++.h>
#include <chrono>
using namespace std;

int d, t;
double eps = 1e-9;
int numColors = 1; // set dynamically after reading colors

struct Point
{
    vector<double> coords; // d + t dimensions combined
    int index;
    int color_id;
};

struct QueryBox
{
    vector<pair<double, double>> bounds;
};

struct FeatureBox
{
    vector<pair<double, double>> bounds;
};

struct RangeTreeNode
{
    int level;
    double split_val;
    int subtree_size;
    vector<int> color_count;
    RangeTreeNode *left;
    RangeTreeNode *right;
    RangeTreeNode *next_tree;
    Point *min_point;
    Point *curr_point;

    RangeTreeNode(int l, double val, int subcnt, RangeTreeNode *le, RangeTreeNode *ri,
                  RangeTreeNode *ne, Point *cp, Point *mp, vector<int> v)
        : level(l), split_val(val), subtree_size(subcnt), left(le), right(ri),
          next_tree(ne), curr_point(cp), min_point(mp), color_count(v) {}
};

vector<Point> points;

struct HashPair
{
    template <class T1, class T2>
    size_t operator()(const std::pair<T1, T2> &p) const
    {
        auto hash1 = std::hash<T1>{}(p.first);
        auto hash2 = std::hash<T2>{}(p.second);
        return hash1 ^ (hash2 + 0x9e3779b9 + (hash1 << 6) + (hash1 >> 2));
    }
};

struct HashVectorOfPairs
{
    size_t operator()(const std::vector<std::pair<double, double>> &v) const
    {
        size_t seed = v.size();
        for (const auto &p : v)
        {
            seed ^= (HashPair{}(p) + 0x9e3779b9 + (seed << 6) + (seed >> 2));
        }
        return seed;
    }
};

static inline bool inRangeBox(const Point &p, const QueryBox &box)
{
    for (int i = 0; i < d; ++i)
    {
        if (p.coords[i] < box.bounds[i].first || p.coords[i] > box.bounds[i].second)
            return false;
    }
    return true;
}

static inline bool inFeatureBox(const Point &p, const FeatureBox &box)
{
    for (int i = 0; i < t; ++i)
    {
        if (p.coords[d + i] < box.bounds[i].first || p.coords[d + i] > box.bounds[i].second)
            return false;
    }
    return true;
}

static inline bool min_all_masked(Point *p1, Point *p2, int mask)
{
    for (int i = d; i >= 0; i--)
    {
        int dim_idx = d + i;
        double val1 = p1->coords[dim_idx];
        double val2 = p2->coords[dim_idx];
        if (((mask >> i) & 1) == 0)
        {
            if (val1 < val2)
                return true;
            if (val1 > val2)
                return false;
        }
        else
        {
            if (val1 > val2)
                return true;
            if (val1 < val2)
                return false;
        }
    }
    return false;
}

static inline bool min_all(Point *p1, Point *p2)
{
    for (int i = t - 1; i >= 0; i--)
    {
        if (p1->coords[d + i] < p2->coords[d + i])
            return true;
        if (p1->coords[d + i] > p2->coords[d + i])
            return false;
    }
    return false;
}

RangeTreeNode *build_tree(vector<Point *> &pts, int level, int last_dim, int total_levels)
{
    if (pts.empty())
        return nullptr;

    int dim = level >= last_dim ? level + 1 : level;
    if (level == total_levels - 1)
        dim = last_dim;

    sort(pts.begin(), pts.end(), [&](Point *a, Point *b)
         { return a->coords[dim] < b->coords[dim]; });

    int mid = (int)pts.size() / 2;
    Point *curr = pts[mid];
    double split_val = curr->coords[dim];

    vector<Point *> left_pts(pts.begin(), pts.begin() + mid);
    vector<Point *> right_pts(pts.begin() + mid + 1, pts.end());

    RangeTreeNode *left = build_tree(left_pts, level, last_dim, total_levels);
    RangeTreeNode *right = build_tree(right_pts, level, last_dim, total_levels);
    RangeTreeNode *next = nullptr;
    if (level < total_levels - 1)
        next = build_tree(pts, level + 1, last_dim, total_levels);

    vector<int> v(numColors, 0);
    if (left)
        for (int i = 0; i < numColors; i++)
            v[i] += left->color_count[i];
    if (right)
        for (int i = 0; i < numColors; i++)
            v[i] += right->color_count[i];
    if (curr->color_id >= 0 && curr->color_id < numColors)
        v[curr->color_id]++;

    Point *min_p = curr;
    if (left && min_all(left->min_point, min_p))
        min_p = left->min_point;
    if (right && min_all(right->min_point, min_p))
        min_p = right->min_point;
    if (next && min_all(next->min_point, min_p))
        min_p = next->min_point;

    return new RangeTreeNode(level, split_val, (int)pts.size(), left, right, next, curr, min_p, v);
}

static inline bool inRangeQ(Point *p, QueryBox &Q, int D, int T)
{
    for (int i = 0; i < D; ++i)
    {
        if (p->coords[i] < Q.bounds[i].first || p->coords[i] > Q.bounds[i].second)
            return false;
    }
    return true;
}

static inline bool inRangeQExpand(Point *p, QueryBox &Q, int D, int T, int last_dim)
{
    for (int i = 0; i < D; ++i)
    {
        if (i == last_dim)
            continue;
        if (p->coords[i] < Q.bounds[i].first || p->coords[i] > Q.bounds[i].second)
            return false;
    }
    return true;
}

static inline bool inRangeF(Point *p, FeatureBox &F, int D, int T)
{
    for (int i = D; i < D + T; ++i)
    {
        if (p->coords[i] < F.bounds[i].first || p->coords[i] > F.bounds[i].second)
            return false;
    }
    return true;
}

RangeTreeNode *findSplitNode(RangeTreeNode *root, double low, double high)
{
    while (root && (root->left || root->right) && (high < root->split_val || low > root->split_val))
    {
        if (high < root->split_val)
            root = root->left;
        else
            root = root->right;
    }
    return root;
}

void dDRangeQuery(RangeTreeNode *root, QueryBox &Q, FeatureBox &F, int level, int total_levels,
                  int last_dim, int d, int t, vector<Point *> &out, int &size, vector<int> &color_count)
{
    if (!root || level >= total_levels)
        return;

    int dim;
    if (level < d)
    {
        dim = level >= last_dim ? level + 1 : level;
        if (level == total_levels - 1)
            dim = last_dim;
    }
    else
    {
        dim = level - d;
    }

    double low = (level < d) ? Q.bounds[dim].first : F.bounds[dim].first;
    double high = (level < d) ? Q.bounds[dim].second : F.bounds[dim].second;

    RangeTreeNode *split = findSplitNode(root, low, high);
    if (!split)
        return;

    if (inRangeQ(split->curr_point, Q, d, t))
    {
        out.push_back(split->curr_point);
        size++;
        if (split->curr_point->color_id >= 0 && split->curr_point->color_id < numColors)
            color_count[split->curr_point->color_id]++;
    }

    // left path
    RangeTreeNode *v = split->left;
    while (v)
    {
        if (low <= v->split_val)
        {
            if (inRangeQ(v->curr_point, Q, d, t))
            {
                out.push_back(v->curr_point);
                size++;
                if (v->curr_point->color_id >= 0 && v->curr_point->color_id < numColors)
                    color_count[v->curr_point->color_id]++;
            }
            if (level == total_levels - 1)
            {
                if (v->right)
                    size += v->right->subtree_size;
                if (v->right)
                    for (int i = 0; i < numColors; i++)
                        color_count[i] += v->right->color_count[i];
            }
            if (v->right && v->right->next_tree)
            {
                dDRangeQuery(v->right->next_tree, Q, F, level + 1, total_levels, last_dim, d, t, out, size, color_count);
            }
            v = v->left;
        }
        else
        {
            v = v->right;
        }
    }

    // right path
    v = split->right;
    while (v)
    {
        if (high >= v->split_val)
        {
            if (inRangeQ(v->curr_point, Q, d, t))
            {
                out.push_back(v->curr_point);
                size++;
                if (v->curr_point->color_id >= 0 && v->curr_point->color_id < numColors)
                    color_count[v->curr_point->color_id]++;
            }
            if (level == total_levels - 1)
            {
                if (v->left)
                    size += v->left->subtree_size;
                if (v->left)
                    for (int i = 0; i < numColors; i++)
                        color_count[i] += v->left->color_count[i];
            }
            if (v->left && v->left->next_tree)
            {
                dDRangeQuery(v->left->next_tree, Q, F, level + 1, total_levels, last_dim, d, t, out, size, color_count);
            }
            v = v->right;
        }
        else
        {
            v = v->left;
        }
    }
}

void side_expand_2(RangeTreeNode *root, QueryBox &Q, FeatureBox &F, int level, int total_levels,
                   int last_dim, int d, int t, vector<Point *> &out, int &size, vector<int> &color_count)
{
    if (!root || level >= total_levels)
        return;

    int dim = (level < d) ? level : (level - d);
    if (level < d)
    {
        dim = level >= last_dim ? level + 1 : level;
        if (level == total_levels - 1)
            dim = last_dim;
    }
    double low = (level < d) ? Q.bounds[dim].first : F.bounds[dim].first;
    double high = (level < d) ? Q.bounds[dim].second : F.bounds[dim].second;

    RangeTreeNode *split = findSplitNode(root, low, high);
    if (!split)
        return;

    if (inRangeQExpand(split->curr_point, Q, d, t, last_dim))
    {
        out.push_back(split->curr_point);
        size++;
        if (split->curr_point->color_id >= 0 && split->curr_point->color_id < numColors)
            color_count[split->curr_point->color_id]++;
    }

    RangeTreeNode *v = split->left;
    while (v)
    {
        if (low <= v->split_val)
        {
            if (inRangeQExpand(v->curr_point, Q, d, t, last_dim))
            {
                out.push_back(v->curr_point);
                size++;
                if (v->curr_point->color_id >= 0 && v->curr_point->color_id < numColors)
                    color_count[v->curr_point->color_id]++;
            }
            if (level == total_levels - 1)
            {
                if (v->right)
                    size += v->right->subtree_size;
                if (v->right)
                    for (int i = 0; i < numColors; i++)
                        color_count[i] += v->right->color_count[i];
            }
            if (v->right && v->right->next_tree)
            {
                side_expand_2(v->right->next_tree, Q, F, level + 1, total_levels, last_dim, d, t, out, size, color_count);
            }
            v = v->left;
        }
        else
        {
            if (level == total_levels - 1)
            {
                out.push_back(v->curr_point);
                size++;
                if (v->curr_point->color_id >= 0 && v->curr_point->color_id < numColors)
                    color_count[v->curr_point->color_id]++;
            }
            v = v->right;
        }
    }

    v = split->right;
    while (v)
    {
        if (high >= v->split_val)
        {
            if (inRangeQExpand(v->curr_point, Q, d, t, last_dim))
            {
                out.push_back(v->curr_point);
                size++;
                if (v->curr_point->color_id >= 0 && v->curr_point->color_id < numColors)
                    color_count[v->curr_point->color_id]++;
            }
            if (level == total_levels - 1)
            {
                if (v->left)
                    size += v->left->subtree_size;
                if (v->left)
                    for (int i = 0; i < numColors; i++)
                        color_count[i] += v->left->color_count[i];
            }
            if (v->left && v->left->next_tree)
            {
                side_expand_2(v->left->next_tree, Q, F, level + 1, total_levels, last_dim, d, t, out, size, color_count);
            }
            v = v->right;
        }
        else
        {
            if (level == total_levels - 1)
            {
                out.push_back(v->curr_point);
                size++;
                if (v->curr_point->color_id >= 0 && v->curr_point->color_id < numColors)
                    color_count[v->curr_point->color_id]++;
            }
            v = v->left;
        }
    }
}

vector<FeatureBox> splitFeatureBoxMasked(const FeatureBox &fbox, const Point &p, int mask)
{
    vector<FeatureBox> subBoxes;
    for (int i = 0; i < d; i++)
    {
        FeatureBox box = fbox;
        for (int j = 0; j < i; j++)
        {
            if (((mask >> j) & 1) == 0)
                box.bounds[j].first = max(p.coords[j], box.bounds[j].first);
            else
                box.bounds[j].second = min(p.coords[j], box.bounds[j].second);
        }
        if (((mask >> i) & 1) == 0)
            box.bounds[i].second = min(p.coords[i] - 1e-9, box.bounds[i].second);
        else
            box.bounds[i].first = max(p.coords[i] + 1e-9, box.bounds[i].first);

        bool valid = true;
        for (int j = 0; j < d; j++)
            if (box.bounds[j].first > box.bounds[j].second)
            {
                valid = false;
                break;
            }
        for (int j = 0; j < d; j++)
            if (box.bounds[j].first == fbox.bounds[j].first && box.bounds[j].second == fbox.bounds[j].second)
            {
                valid = false;
                break;
            }
        if (valid)
            subBoxes.push_back(box);
    }
    return subBoxes;
}

vector<FeatureBox> splitFeatureBox(const FeatureBox &fbox, const Point &p)
{
    vector<FeatureBox> subBoxes;
    for (int i = 0; i < t - 1; i++)
    {
        FeatureBox box = fbox;
        for (int j = 0; j < i; j++)
            box.bounds[j].first = max(p.coords[d + j], box.bounds[j].first);
        box.bounds[i].second = min(p.coords[d + i] - 1e-9, box.bounds[i].second);
        box.bounds[t - 1].first = max(p.coords[d + t - 1], box.bounds[t - 1].first);

        bool valid = true;
        for (int j = 0; j < t; j++)
            if (box.bounds[j].first > box.bounds[j].second)
            {
                valid = false;
                break;
            }
        if (valid)
            subBoxes.push_back(box);
    }
    return subBoxes;
}

vector<Point *> findSkyline(RangeTreeNode *root, const QueryBox &qbox, int mask)
{
    FeatureBox fullBox;
    fullBox.bounds.resize(d);
    for (int i = 0; i < d; i++)
        fullBox.bounds[i] = {qbox.bounds[i].first, qbox.bounds[i].second};
    stack<FeatureBox> Z;
    Z.push(fullBox);

    unordered_set<int> reported;
    vector<Point *> skyline;

    while (!Z.empty())
    {
        FeatureBox R = Z.top();
        Z.pop();
        QueryBox Rq;
        for (int i = 0; i < d; i++)
            Rq.bounds.push_back({R.bounds[i].first, R.bounds[i].second});

        vector<Point *> canonicals;
        int size = 0;
        vector<int> temp(numColors, 0);
        dDRangeQuery(root, const_cast<QueryBox &>(Rq), R, 0, d, d - 1, d, t, canonicals, size, temp);

        Point *best = nullptr;
        for (auto pt : canonicals)
        {
            if (!inRangeBox(*pt, qbox))
                continue;
            if (!best || min_all_masked(pt, best, mask))
                best = pt;
        }

        if (best)
        {
            auto subBoxes = splitFeatureBoxMasked(R, *best, mask);
            for (auto it = subBoxes.rbegin(); it != subBoxes.rend(); ++it)
                Z.push(*it);
            if (reported.count(best->index) == 0)
            {
                skyline.push_back(best);
                reported.insert(best->index);
            }
        }
    }
    return skyline;
}

static inline bool dominates(const Point &p, const Point &q)
{
    bool strictlyBetter = false;
    for (int i = t - 1; i >= 0; i--)
    {
        if (p.coords[d + i] > q.coords[d + i])
            return false;
        if (p.coords[d + i] < q.coords[d + i])
            strictlyBetter = true;
    }
    return strictlyBetter;
}

static inline bool isSamebound(QueryBox &Q1, QueryBox &Q2, int dim)
{
    for (int i = 0; i < dim; i++)
    {
        if (Q1.bounds[i].first != Q2.bounds[i].first)
            return false;
        if (Q1.bounds[i].second != Q2.bounds[i].second)
            return false;
    }
    return true;
}

vector<QueryBox> generateCornerSkylineRanges(const QueryBox &qbox,
                                             const vector<double> &best_lhs,
                                             const vector<double> &best_rhs)
{
    int d = (int)qbox.bounds.size();
    vector<QueryBox> corner_ranges;
    QueryBox copybox = qbox;
    int total_masks = 1 << d;

    for (int mask = 0; mask < total_masks; ++mask)
    {
        QueryBox newbox = qbox;
        for (int i = 0; i < d; ++i)
        {
            if ((mask >> i) & 1)
            {
                newbox.bounds[i].first = newbox.bounds[i].second + 1e-9;
                newbox.bounds[i].second = best_rhs[i];
            }
            else
            {
                newbox.bounds[i].second = newbox.bounds[i].first - 1e-9;
                newbox.bounds[i].first = best_lhs[i];
            }
        }
        bool valid = true;
        for (int i = 0; i < d; ++i)
        {
            if (newbox.bounds[i].first > newbox.bounds[i].second)
            {
                valid = false;
                break;
            }
        }
        if (isSamebound(newbox, copybox, d))
            valid = false;
        if (!valid)
            continue;
        corner_ranges.push_back(newbox);
    }
    return corner_ranges;
}

QueryBox intersectBoxes(const QueryBox &Q, const QueryBox &neighbor, int d, bool &valid)
{
    QueryBox intersectBox;
    intersectBox.bounds.resize(d);
    valid = true;
    for (int i = 0; i < d; i++)
    {
        double low = std::max(Q.bounds[i].first, neighbor.bounds[i].first);
        double high = std::min(Q.bounds[i].second, neighbor.bounds[i].second);
        if (low > high)
        {
            valid = false;
            break;
        }
        intersectBox.bounds[i] = {low, high};
    }
    return intersectBox;
}

QueryBox unionBoxes(const QueryBox &Q, const QueryBox &neighbor, int d)
{
    QueryBox unionBox;
    unionBox.bounds.resize(d);
    for (int i = 0; i < d; i++)
    {
        double low = std::min(Q.bounds[i].first, neighbor.bounds[i].first);
        double high = std::max(Q.bounds[i].second, neighbor.bounds[i].second);
        unionBox.bounds[i] = {low, high};
    }
    return unionBox;
}

// ---- Fairness helper for K colors ----
static inline bool isFairK(const vector<int> &counts)
{
    int ref = -1;
    for (int c = 1; c < (int)counts.size(); ++c)
    {
        if (ref == -1)
            ref = counts[c];
        else if (counts[c] != ref)
            return false;
    }
    return ref > 0; // non-empty
}

// BFS over boxes; returns [d pairs for range] + [1 pair with sim]
// NOTE: now takes an explicit initial box to start from.
vector<pair<double, double>> bfs(vector<RangeTreeNode *> &trees,
                                 QueryBox &qbox,
                                 RangeTreeNode *root,
                                 const QueryBox &init_box)
{
    FeatureBox fullBox;
    unordered_set<std::vector<std::pair<double, double>>, HashVectorOfPairs> boxSet;

    fullBox.bounds.resize(t);
    for (int i = 0; i < t; i++)
    {
        fullBox.bounds[i] = {-numeric_limits<double>::infinity(), numeric_limits<double>::infinity()};
    }

    priority_queue<vector<double>> pq; // {sim, pos}
    vector<QueryBox> boxes;

    // Start from the user-provided initial box (already clamped to domain in main)
    boxes.push_back(init_box);
    boxSet.insert(init_box.bounds);
    pq.push({1.0, 0.0}); // initial sim = 1.0, index = 0

    // Base color counts inside original domain qbox (unchanged)
    int base_size = 0;
    vector<Point *> base_canonicals;
    vector<int> base_color_count(numColors, 0);
    dDRangeQuery(root, const_cast<QueryBox &>(qbox), fullBox, 0, d, d - 1, d, t,
                 base_canonicals, base_size, base_color_count);

    // cout << "Base Color count:\n"; for (int i = 1; i < numColors; i++) cout << base_color_count[i] << " "; cout << "\n";

    QueryBox intersect_box, union_box;

    while (!pq.empty())
    {
        auto priority = pq.top();
        pq.pop();
        double sim = priority[0];
        int pos = (int)priority[1];

        QueryBox box = boxes[pos];

        // Current counts in this box
        vector<Point *> cur_canonicals;
        int subtree_size = 0;
        vector<int> color_count(numColors, 0);
        dDRangeQuery(root, const_cast<QueryBox &>(box), fullBox, 0, d, d - 1, d, t,
                     cur_canonicals, subtree_size, color_count);

        // Exact recount for fairness check
        vector<int> exact(numColors, 0);
        for (int i = 0; i < (int)points.size(); i++)
        {
            if (inRangeBox(points[i], box))
            {
                int c = points[i].color_id;
                if (c >= 0 && c < numColors)
                    exact[c]++;
            }
        }

        if (isFairK(exact))
        {
            vector<pair<double, double>> res;
            for (int i = 0; i < d; i++)
                res.push_back({box.bounds.at(i).first, box.bounds.at(i).second});
            res.push_back({sim, sim});
            return res;
        }

        vector<double> lhs_expands, rhs_expands;

        for (int dim = 0; dim < d; ++dim)
        {
            vector<Point *> out;
            int size = 0;
            vector<int> temp(numColors, 0);
            side_expand_2(trees[dim], box, fullBox, 0, d, dim, d, t, out, size, temp);

            Point *lhs_shrink_pt_ptr = nullptr;
            Point *rhs_shrink_pt_ptr = nullptr;
            Point *lhs_expand_pt_ptr = nullptr;
            Point *rhs_expand_pt_ptr = nullptr;

            double lhs_bound = -1e16, rhs_bound = 1e16;
            double lhs_shrink_pt = box.bounds[dim].second, rhs_shrink_pt = box.bounds[dim].first;

            for (auto pt : out)
            {
                if (pt->coords[dim] < box.bounds[dim].first)
                {
                    if (pt->coords[dim] > lhs_bound)
                    {
                        lhs_bound = pt->coords[dim];
                        lhs_expand_pt_ptr = pt;
                    }
                }
                else if (pt->coords[dim] > box.bounds[dim].second)
                {
                    if (pt->coords[dim] < rhs_bound)
                    {
                        rhs_bound = pt->coords[dim];
                        rhs_expand_pt_ptr = pt;
                    }
                }
                else
                {
                    if (pt->coords[dim] < lhs_shrink_pt)
                    {
                        lhs_shrink_pt = pt->coords[dim];
                        lhs_shrink_pt_ptr = pt;
                    }
                    if (pt->coords[dim] > rhs_shrink_pt)
                    {
                        rhs_shrink_pt = pt->coords[dim];
                        rhs_shrink_pt_ptr = pt;
                    }
                }
            }

            QueryBox lhs_expanded_box = box, rhs_expanded_box = box, lhs_shrink_box = box, rhs_shrink_box = box;

            lhs_expands.push_back(lhs_bound);
            rhs_expands.push_back(rhs_bound);

            // LHS expand
            lhs_expanded_box.bounds[dim].first = lhs_bound;
            bool valid1 = false;
            intersect_box = intersectBoxes(qbox, lhs_expanded_box, d, valid1);
            union_box = unionBoxes(qbox, lhs_expanded_box, d);
            if (!isSamebound(lhs_expanded_box, box, d) && valid1)
            {
                int union_size = 0, intersection_size = 0;
                vector<Point *> tmp;
                vector<int> ic(numColors, 0), uc(numColors, 0);
                dDRangeQuery(root, const_cast<QueryBox &>(intersect_box), fullBox, 0, d, d - 1, d, t, tmp, intersection_size, ic);
                dDRangeQuery(root, const_cast<QueryBox &>(union_box), fullBox, 0, d, d - 1, d, t, tmp, union_size, uc);
                double cor_sim = (double)(intersection_size) / (double)(union_size);
                if (cor_sim < sim && boxSet.count(lhs_expanded_box.bounds) == 0)
                {
                    boxes.push_back(lhs_expanded_box);
                    boxSet.insert(lhs_expanded_box.bounds);
                    pq.push({cor_sim, (double)(boxes.size() - 1)});
                }
            }

            // RHS expand
            rhs_expanded_box.bounds[dim].second = rhs_bound;
            bool valid2 = false;
            intersect_box = intersectBoxes(qbox, rhs_expanded_box, d, valid2);
            union_box = unionBoxes(qbox, rhs_expanded_box, d);
            if (!isSamebound(rhs_expanded_box, box, d) && valid2)
            {
                int union_size = 0, intersection_size = 0;
                vector<Point *> tmp;
                vector<int> ic(numColors, 0), uc(numColors, 0);
                dDRangeQuery(root, const_cast<QueryBox &>(intersect_box), fullBox, 0, d, d - 1, d, t, tmp, intersection_size, ic);
                dDRangeQuery(root, const_cast<QueryBox &>(union_box), fullBox, 0, d, d - 1, d, t, tmp, union_size, uc);
                double cor_sim = (double)(intersection_size) / (double)(union_size);
                if (cor_sim < sim && boxSet.count(rhs_expanded_box.bounds) == 0)
                {
                    boxes.push_back(rhs_expanded_box);
                    boxSet.insert(rhs_expanded_box.bounds);
                    pq.push({cor_sim, (double)(boxes.size() - 1)});
                }
            }

            // LHS shrink
            lhs_shrink_box.bounds[dim].first = lhs_shrink_pt;
            bool valid3 = false;
            intersect_box = intersectBoxes(qbox, lhs_shrink_box, d, valid3);
            union_box = unionBoxes(qbox, lhs_shrink_box, d);
            if (!isSamebound(lhs_shrink_box, box, d) && valid3)
            {
                int union_size = 0, intersection_size = 0;
                vector<Point *> tmp;
                vector<int> ic(numColors, 0), uc(numColors, 0);
                dDRangeQuery(root, const_cast<QueryBox &>(intersect_box), fullBox, 0, d, d - 1, d, t, tmp, intersection_size, ic);
                dDRangeQuery(root, const_cast<QueryBox &>(union_box), fullBox, 0, d, d - 1, d, t, tmp, union_size, uc);
                union_size += base_size - intersection_size;
                double cor_sim = (double)(intersection_size) / (double)(union_size);
                if (cor_sim < sim && boxSet.count(lhs_shrink_box.bounds) == 0)
                {
                    boxes.push_back(lhs_shrink_box);
                    boxSet.insert(lhs_shrink_box.bounds);
                    pq.push({cor_sim, (double)(boxes.size() - 1)});
                }
            }

            // RHS shrink
            rhs_shrink_box.bounds[dim].second = rhs_shrink_pt;
            bool valid4 = false;
            intersect_box = intersectBoxes(qbox, rhs_shrink_box, d, valid4);
            union_box = unionBoxes(qbox, rhs_shrink_box, d);
            if (!isSamebound(rhs_shrink_box, box, d) && valid4)
            {
                int union_size = 0, intersection_size = 0;
                vector<Point *> tmp;
                vector<int> ic(numColors, 0), uc(numColors, 0);
                dDRangeQuery(root, const_cast<QueryBox &>(intersect_box), fullBox, 0, d, d - 1, d, t, tmp, intersection_size, ic);
                dDRangeQuery(root, const_cast<QueryBox &>(union_box), fullBox, 0, d, d - 1, d, t, tmp, union_size, uc);
                union_size += base_size - intersection_size;
                double cor_sim = (double)(intersection_size) / (double)(union_size);
                if (cor_sim < sim && boxSet.count(rhs_shrink_box.bounds) == 0)
                {
                    boxes.push_back(rhs_shrink_box);
                    boxSet.insert(rhs_shrink_box.bounds);
                    pq.push({cor_sim, (double)(boxes.size() - 1)});
                }
            }
        }

        // Corner expansions via skyline
        vector<QueryBox> corner_boxes = generateCornerSkylineRanges(box, lhs_expands, rhs_expands);
        for (int i = 0; i < (int)corner_boxes.size(); i++)
        {
            const auto &corner_box = corner_boxes[i];
            vector<Point *> skyline = findSkyline(root, corner_box, i);
            for (Point *pt : skyline)
            {
                QueryBox neighbor = box;
                for (int j = 0; j < d; ++j)
                {
                    neighbor.bounds[j].first = min(neighbor.bounds[j].first, pt->coords[j]);
                    neighbor.bounds[j].second = max(neighbor.bounds[j].second, pt->coords[j]);
                }
                bool valid = false;
                intersect_box = intersectBoxes(qbox, neighbor, d, valid);
                union_box = unionBoxes(qbox, neighbor, d);
                if (!isSamebound(neighbor, box, d) && valid)
                {
                    int union_size = 0, intersection_size = 0;
                    vector<Point *> tmp;
                    vector<int> ic(numColors, 0), uc(numColors, 0);
                    dDRangeQuery(root, const_cast<QueryBox &>(intersect_box), fullBox, 0, d, d - 1, d, t, tmp, intersection_size, ic);
                    dDRangeQuery(root, const_cast<QueryBox &>(union_box), fullBox, 0, d, d - 1, d, t, tmp, union_size, uc);
                    union_size += base_size - intersection_size;
                    double cor_sim = (double)(intersection_size) / (double)(union_size);
                    if (cor_sim < sim && boxSet.count(neighbor.bounds) == 0)
                    {
                        boxes.push_back(neighbor);
                        boxSet.insert(neighbor.bounds);
                        pq.push({cor_sim, (double)(boxes.size() - 1)});
                    }
                }
            }
        }
    }

    return {}; // no fair box found
}

// ------------------- main -------------------
int main(int argc, char **argv)
{
    using std::ifstream;
    using std::istream;
    using std::string;

    // read from argv[1] (or "-" for stdin), else fallback to "points.txt"
    string fname = "points.txt";
    if (argc >= 2)
        fname = argv[1];
    else
    {
        std::cerr << "No input file arg given. Falling back to 'points.txt'.\n";
        std::cerr << "Usage: " << (argv[0] ? argv[0] : "prog") << " <input-file>|-\n";
    }

    istream *pin = nullptr;
    ifstream fin;
    if (fname == "-")
        pin = &std::cin;
    else
    {
        fin.open(fname);
        if (!fin)
        {
            std::cerr << "Error: cannot open file '" << fname << "'.\n";
            return 1;
        }
        pin = &fin;
    }
    istream &IN = *pin;

    // input
    int n;
    IN >> n;
    IN >> d >> t;

    points.clear();
    points.resize(n);
    for (int i = 0; i < n; i++)
    {
        points[i].coords.resize(d + t);
        for (int j = 0; j < d + t; j++)
            IN >> points[i].coords[j];
        points[i].index = i;
    }

    // labels, infer numColors
    int maxColorId = 0;
    for (int i = 0; i < n; i++)
    {
        IN >> points[i].color_id;
        maxColorId = max(maxColorId, points[i].color_id);
    }
    numColors = std::max(2, maxColorId + 1);

    // build (timed)
    using clock_t = std::chrono::steady_clock;
    auto t_build_start = clock_t::now();

    std::vector<Point *> point_ptrs(n);
    for (int i = 0; i < n; ++i)
        point_ptrs[i] = &points[i];

    RangeTreeNode *root = build_tree(point_ptrs, 0, d + t - 1, d + t);
    std::vector<RangeTreeNode *> trees(d);
    for (int i = 0; i < d; i++)
        trees[i] = build_tree(point_ptrs, 0, i, d);

    auto t_build_end = clock_t::now();
    auto build_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_build_end - t_build_start).count();

    // domain (d lines)
    QueryBox qbox;
    qbox.bounds.resize(d);
    for (int i = 0; i < d; i++)
        IN >> qbox.bounds[i].first >> qbox.bounds[i].second;

    // initial box from user (d lines), then clamp to domain
    QueryBox init_box;
    init_box.bounds.resize(d);
    for (int i = 0; i < d; ++i)
    {
        IN >> init_box.bounds[i].first >> init_box.bounds[i].second;
        if (init_box.bounds[i].first > init_box.bounds[i].second)
            std::swap(init_box.bounds[i].first, init_box.bounds[i].second);
        init_box.bounds[i].first = std::max(init_box.bounds[i].first, qbox.bounds[i].first);
        init_box.bounds[i].second = std::min(init_box.bounds[i].second, qbox.bounds[i].second);
    }

    // query (timed) starting from user-provided seed
    auto t_query_start = clock_t::now();
    std::vector<std::pair<double, double>> ans = bfs(trees, qbox, root, init_box);
    auto t_query_end = clock_t::now();
    auto query_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_query_end - t_query_start).count();

    // timings
    std::cout << "Build time (ms): " << build_ms << "\n";
    std::cout << "Query time (ms): " << query_ms << "\n";
    std::cout << "Total time (ms): " << (build_ms + query_ms) << "\n";

    // final summary
    std::cout << "Ans size " << ans.size() << "\n";
    if ((int)ans.size() >= d + 1)
    {
        for (int i = 0; i < d; i++)
            std::cout << "[" << ans[i].first << "," << ans[i].second << "], ";
        std::cout << "\nMax sim: " << ans[d].first << "\n";

        QueryBox final_box;
        final_box.bounds.resize(d);
        for (int i = 0; i < d; i++)
        {
            final_box.bounds[i].first = ans[i].first;
            final_box.bounds[i].second = ans[i].second;
        }

        std::cout << "Points in final box:\n";
        std::vector<int> check(numColors, 0);
        for (int i = 0; i < n; i++)
        {
            if (inRangeQ(&points[i], final_box, d, t))
            {
                int c = points[i].color_id;
                if (c >= 0 && c < numColors)
                    check[c]++;
            }
        }
        std::cout << "Color counts: ";
        for (int c = 1; c < numColors; ++c)
            std::cout << check[c] << " ";
        std::cout << "\n";
    }
    else
    {
        std::cout << "No fair box found.\n";
    }
    return 0;
}
