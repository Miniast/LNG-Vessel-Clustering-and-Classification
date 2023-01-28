#include "Cluster.h"

inline ll get_key(double x, double y)
{
    ll a = (ll)round(10000 * x) + 1600000, b = (ll)round(10000 * y) + 900000;
    return a * 10000000 + b;
}
inline double distance(nodes &x, cluster &y)
{
    return sqrt(pow((x.latitude - y.avg_lati), 2) + pow(x.longtitude - y.avg_longti, 2));
}
inline double distance(nodes &x, nodes &y)
{
    return fabs(x.latitude - y.latitude) + fabs(x.longtitude - y.longtitude);
}
inline double distance(cluster &x, cluster &y)
{
    return sqrt(pow((x.avg_lati - y.avg_lati), 2) + pow(x.avg_longti - y.avg_longti, 2));
}
inline bool cmp(int x, int y)
{
    if (vessel[x].mmsi == vessel[y].mmsi)
        return vessel[x].unix_time < vessel[y].unix_time;
    return vessel[x].mmsi < vessel[y].mmsi;
}
inline bool cmp2(nodes &x, nodes &y)
{
    if (x.longtitude == y.longtitude)
        return x.latitude < y.latitude;
    return x.longtitude < y.longtitude;
}
inline bool cmp3(cluster x, cluster y)
{
    if (x.avg_longti == y.avg_longti)
        return x.avg_lati < y.avg_lati;
    return x.avg_longti < y.avg_longti;
}
void pre_treat()
{
    R int cnt = 0, pre;
    for (int i = 1; i <= m; ++i)
        while (batch[i].number <= 5)
        {
            swap(batch[i], batch[m]);
            m--;
        }

    sort(batch + 1, batch + m + 1, cmp2);

    pre = 1;
    batch[m + 1].longtitude = 1000.0;
    for (int i = 1; i <= m; ++i)
    {
        if (distance(batch[i], batch[i + 1]) < 5e-4)
        {
            batch[i + 1].longtitude = (batch[i].longtitude * batch[i].number + batch[i + 1].longtitude * batch[i + 1].number) / (double)(batch[i + 1].number + batch[i].number);
            batch[i + 1].latitude = (batch[i].latitude * batch[i].number + batch[i + 1].latitude * batch[i + 1].number) / (double)(batch[i + 1].number + batch[i].number);
            batch[i + 1].number += batch[i].number;
            batch[i].number = 0;
        }
        else
        {
            for (int j = pre; j < i; ++j)
                batch[i].node_list.insert(batch[i].node_list.end(), batch[j].node_list.begin(), batch[j].node_list.end());
            pre = i + 1;
        }
    }

    cnt = 0;
    for (int i = 1; i <= m; ++i)
        if (batch[i].number > 0)
            batch[++cnt] = batch[i];
    m = cnt;

    cnt = 0;
    for (int i = 1; i <= m; ++i)
    {
        if (batch[i].number >= 500)
            ord[++cnt] = i;
    }
    K = cnt;

    for (int i = 1; i <= K; ++i)
    {
        clst[i].avg_longti = batch[ord[i]].longtitude;
        clst[i].avg_lati = batch[ord[i]].latitude;
        clst[i].number = 0;
    }
}
void Mini_Batch_K_Means()
{
    R bool find = false;
    R int tar = 0, cntk;
    R double dist, tmp, avg_x, avg_y;
    while (!find)
    {
        times++;
        if (times > 50)
            break;
        for (int i = 1; i <= m; ++i)
        {
            dist = inf;
            for (int j = 1; j <= K; ++j)
            {
                tmp = distance(batch[i], clst[j]);
                if (dist > tmp)
                {
                    tar = j;
                    dist = tmp;
                }
            }
            clst[tar].batch_list.push_back(i);
            clst[tar].number += batch[i].number;
            clst[tar].sum_longti += batch[i].longtitude * batch[i].number;
            clst[tar].sum_lati += batch[i].latitude * batch[i].number;
        }
        find = true;
        for (int i = 1; i <= K; ++i)
            if (clst[i].number > 0)
            {
                avg_x = clst[i].sum_longti / clst[i].number;
                avg_y = clst[i].sum_lati / clst[i].number;
                if (fabs(avg_x - clst[i].avg_longti) > 1e-6 || fabs(avg_y - clst[i].avg_lati) > 1e-6)
                    find = false;
                clst[i].avg_longti = avg_x;
                clst[i].avg_lati = avg_y;
            }
            else
                clst[i].avg_longti = clst[i].avg_lati = inf;
        if (!find || times == 1)
            for (int i = 1; i <= K; ++i)
            {
                clst[i].sum_longti = clst[i].sum_lati = 0.0;
                clst[i].number = 0;
                clst[i].batch_list.clear();
            }
    }

    sort(clst + 1, clst + K + 1, cmp3);

    for (int i = 1; i <= K; ++i)
    {
        if (distance(clst[i], clst[i + 1]) < 5e-3)
        {
            clst[i + 1].avg_longti = (clst[i].avg_longti * clst[i].number + clst[i + 1].avg_longti * clst[i + 1].number) / (double)(clst[i].number + clst[i + 1].number);
            clst[i + 1].avg_lati = (clst[i].avg_lati * clst[i].number + clst[i + 1].avg_lati * clst[i + 1].number) / (double)(clst[i].number + clst[i + 1].number);
            clst[i + 1].number += clst[i].number;
            clst[i + 1].batch_list.insert(clst[i + 1].batch_list.end(), clst[i].batch_list.begin(), clst[i].batch_list.end());
            clst[i].number = 0;
        }
    }
    cntk = 0;
    for (int i = 1; i <= K; ++i)
        if (clst[i].number > 0)
            clst[++cntk] = clst[i];
    K = cntk;
}
void classfication()
{

    int lst, cnt1, cnt2, cntk;
    for (int i = 1; i <= K; ++i)
    {
        lst = 0;
        cnt1 = 0;
        cnt2 = 0;

        int lb = clst[i].batch_list.size();
        for (int j = 0; j < lb; ++j)
        {
            int cur = clst[i].batch_list[j], ln = batch[cur].node_list.size();
            for (int k = 0; k < ln; ++k)
            {
                int now = batch[cur].node_list[k];
                if (vessel[now].draft > 50 && vessel[now].draft < 300)
                    mylist[++lst] = now;
            }
        }
        sort(mylist + 1, mylist + lst + 1, cmp);

        int d_type = 0;
        for (int j = 2; j <= lst; ++j)
        {
            if (vessel[mylist[j]].mmsi == vessel[mylist[j - 1]].mmsi && vessel[mylist[j]].unix_time - vessel[mylist[j - 1]].unix_time < 10800)
            {
                if (vessel[mylist[j]].draft > vessel[mylist[j - 1]].draft)
                {
                    cnt2 += vessel[mylist[j]].draft - vessel[mylist[j - 1]].draft;
                    if (cnt2 >= 10)
                    {
                        d_type = 2;
                        break;
                    }
                }
                else if (vessel[mylist[j]].draft < vessel[mylist[j - 1]].draft)
                {
                    cnt1 += vessel[mylist[j - 1]].draft - vessel[mylist[j]].draft;
                    if (cnt1 >= 10)
                        d_type = 1;
                }

                else
                    d_type = max(d_type, 0);
            }
        }
        clst[i].type = d_type;
    }
    cntk = 0;
    for (int i = 1; i <= K; ++i)
        if ((clst[i].number >= 500 && clst[i].type == 0) || (clst[i].number >= 2000 && clst[i].type > 0))
            clst[++cntk] = clst[i];
    K = cntk;
}
int main()
{
    clock_t start_time = clock();
    FILE *fp = fopen("lng_results_list.json", "w");

    R int mmsi, unix_time, status, velocity, draft, pos;
    R ll key;
    R double u, v;

    io::CSVReader<7> in("lng2.csv");
    while (in.read_row(mmsi, unix_time, status, velocity, u, v, draft))
    {
        vessel[++n].mmsi = mmsi;
        vessel[n].unix_time = unix_time;
        vessel[n].draft = draft;
        if (!(status == 1 || status == 5 || status == 15) || velocity != 0)
            n--;
        else
        {
            key = get_key(u, v);
            if (mp.find(key) == mp.end())
            {
                pos = ++m;
                mp[key] = pos;
                batch[pos].longtitude = round(u * 10000.0) / 10000.0;
                batch[pos].latitude = round(v * 10000.0) / 10000.0;
            }
            else
                pos = mp[key];
            batch[pos].node_list.push_back(n);
            batch[pos].number++;
        }
    }

    pre_treat();

    Mini_Batch_K_Means();

    classfication();

    int insta = 0, outsta = 0;
    fprintf(fp, "[\n");
    for (int i = 1; i <= K; ++i)
    {
        fprintf(fp, "    {");
        fprintf(fp, "\"code\":%d,\"latitude\":%lf,\"longtitude\":%lf,", i, clst[i].avg_lati, clst[i].avg_longti);
        if (clst[i].type > 0)
        {
            fprintf(fp, "\"isLNG\":true,");
            if (clst[i].type == 1)
                fprintf(fp, "\"IN\":true}");
            else
                fprintf(fp, "\"IN\":false}");
        }
        else
            fprintf(fp, "\"isLNG\":true,\"IN\":None}");
        if (i < K)
            fprintf(fp, ",");
        fprintf(fp, "\n");
        if (clst[i].type == 1)
            insta++;
        else if (clst[i].type == 2)
            outsta++;
    }
    fprintf(fp, "]\n");
    fclose(fp);
    double runtime = (double)(clock() - start_time) / CLOCKS_PER_SEC;
    printf("Runtime = %.2lfs\nTerminal Number = %d, Total Number = %d\n\n", runtime, insta + outsta, K);
    printf("Receiving Terminal:%d\n", insta);
    printf("Export Terminal:%d\n", outsta);
    printf("Anchorage ground:%d\n\n", K - insta - outsta);
    printf("Terminal information are in lng_results_list.json\n");
    return 0;
}

