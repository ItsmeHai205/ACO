

#include <bits/stdc++.h>
using namespace std;

struct City { double x, y; };

double euclid(const City &a, const City &b) {
    double dx = a.x - b.x, dy = a.y - b.y;
    return sqrt(dx*dx + dy*dy);
}

// 2-opt local improvement
double tour_length(const vector<int>& tour, const vector<vector<double>>& dist) {
    double s = 0;
    int n = tour.size();
    for (int i = 0; i < n; ++i) {
        s += dist[tour[i]][tour[(i+1)%n]];
    }
    return s;
}

bool two_opt_once(vector<int>& tour, const vector<vector<double>>& dist) {
    int n = tour.size();
    for (int i = 0; i < n-1; ++i) {
        for (int k = i+1; k < n; ++k) {
            int a = tour[i];
            int b = tour[(i+1)%n];
            int c = tour[k];
            int d = tour[(k+1)%n];
            double delta = dist[a][c] + dist[b][d] - dist[a][b] - dist[c][d];
            if (delta < -1e-9) {
                reverse(tour.begin() + i + 1, tour.begin() + k + 1);
                return true;
            }
        }
    }
    return false;
}

void two_opt(vector<int>& tour, const vector<vector<double>>& dist) {
    while (two_opt_once(tour, dist)) {}
}

int roulette_choice(const vector<double>& probs, mt19937 &rng) {
    double sum = 0;
    for (double v: probs) sum += v;
    if (sum <= 0) {
        vector<int> idx;
        for (int i = 0; i < (int)probs.size(); ++i) if (probs[i] > 0) idx.push_back(i);
        if (idx.empty()) return rng() % probs.size();
        return idx[rng() % idx.size()];
    }
    uniform_real_distribution<double> unif(0.0, sum);
    double r = unif(rng);
    double cum = 0;
    for (int i = 0; i < (int)probs.size(); ++i) {
        cum += probs[i];
        if (r <= cum) return i;
    }
    return probs.size()-1;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) {
        cerr << "Không đọc được n\n";
        return 0;
    }
    vector<City> cities(n);
    for (int i = 0; i < n; ++i) cin >> cities[i].x >> cities[i].y;


    int m = max(10, n);         // số kiến
    int max_iter = 500;         // số vòng lặp
    double alpha = 1.0;         // ảnh hưởng pheromone
    double beta  = 5.0;         // ảnh hưởng heuristic (1/d)
    double rho   = 0.5;         // tốc độ bay hơi
    double Q     = 100.0;       // hằng số pheromone
    bool use_2opt = true;       // bật 2-opt


    // Tính ma trận khoảng cách và heuristic
    vector<vector<double>> dist(n, vector<double>(n, 0.0));
    vector<vector<double>> eta(n, vector<double>(n, 0.0)); 
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) dist[i][j] = 0;
            else {
                dist[i][j] = euclid(cities[i], cities[j]);
            }
            if (dist[i][j] > 0) eta[i][j] = 1.0 / dist[i][j];
            else eta[i][j] = 0.0;
        }
    }

    // Khởi tạo pheromone 
    auto nearest_neighbor_length = [&]() {
        double best = 1e300;
        for (int start = 0; start < min(n, 5); ++start) {
            vector<char> used(n, 0);
            int cur = start;
            used[cur] = 1;
            double L = 0;
            for (int step = 1; step < n; ++step) {
                int nxt = -1;
                double mind = 1e300;
                for (int j = 0; j < n; ++j) if (!used[j] && dist[cur][j] < mind) {
                    mind = dist[cur][j]; nxt = j;
                }
                if (nxt == -1) break;
                L += mind;
                cur = nxt; used[cur] = 1;
            }
     
            L += dist[cur][start];
            best = min(best, L);
        }
        return best;
    };

    double Lnn = nearest_neighbor_length();
    double tau0 = (n>0) ? (1.0 / (n * Lnn)) : 1e-6;

    vector<vector<double>> tau(n, vector<double>(n, tau0)); 


    random_device rd;
    mt19937 rng(rd());

 
    vector<int> global_best_tour;
    double global_best_len = 1e300;


    for (int iter = 0; iter < max_iter; ++iter) {

        vector<vector<int>> ant_tours(m, vector<int>());
        vector<double> ant_lengths(m, 0.0);

        for (int k = 0; k < m; ++k) {
            vector<int> tour;
            tour.reserve(n);
            vector<char> visited(n, 0);

  
            int start = rng() % n;
            tour.push_back(start);
            visited[start] = 1;
            int cur = start;

 
            for (int step = 1; step < n; ++step) {
          
                vector<double> probs(n, 0.0);
                for (int j = 0; j < n; ++j) if (!visited[j]) {
                    probs[j] = pow(tau[cur][j], alpha) * pow(eta[cur][j], beta);
                }
                int nxt = roulette_choice(probs, rng);
         
                if (visited[nxt]) {
                    vector<int> cand;
                    for (int j = 0; j < n; ++j) if (!visited[j]) cand.push_back(j);
                    nxt = cand[rng() % cand.size()];
                }
                tour.push_back(nxt);
                visited[nxt] = 1;
                cur = nxt;
            }
    
            if (use_2opt) two_opt(tour, dist);
            double L = tour_length(tour, dist);
            ant_tours[k] = tour;
            ant_lengths[k] = L;


            if (L < global_best_len) {
                global_best_len = L;
                global_best_tour = tour;
            }
        }


        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                tau[i][j] *= (1.0 - rho);


        for (int k = 0; k < m; ++k) {
            double deposit = Q / ant_lengths[k];
            const vector<int>& tour = ant_tours[k];
            for (int i = 0; i < n; ++i) {
                int a = tour[i];
                int b = tour[(i+1)%n];
                tau[a][b] += deposit;
                tau[b][a] += deposit; 
            }
        }


        double elite_deposit = Q / global_best_len;
        for (int i = 0; i < n; ++i) {
            int a = global_best_tour[i];
            int b = global_best_tour[(i+1)%n];
            tau[a][b] += elite_deposit;
            tau[b][a] += elite_deposit;
        }

      
        double tau_min = 1e-12, tau_max = 1e12;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                tau[i][j] = min(tau_max, max(tau_min, tau[i][j]));


        if (iter % 50 == 0 || iter == max_iter-1) {
            cerr << "Iter " << iter << " best = " << fixed << setprecision(6) << global_best_len << "\n";
        }
    }


    cout << fixed << setprecision(6);
    cout << "Best length: " << global_best_len << "\n";
    cout << "Tour (0-indexed, n=" << global_best_tour.size() << "):\n"s;
    for (int v : global_best_tour) cout << v << " ";
    cout << "\n";
    return 0;
}
