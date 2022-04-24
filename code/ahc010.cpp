// #define SUBMISSION

// 提出時は標準入力から読み込み。練習中はファイルから読み込み
#ifdef SUBMISSION
#define	INPUT(a) cin >> (a);
#define	OUTPUT(a) cout << (a) << endl;
#define	VISUALIZE(a)
#define Print(x)
#define PrintVec(v)
#else
#define	INPUT(a) fin >> (a);
#define	OUTPUT(a) fout << (a) << endl;
#define	VISUALIZE(a) cout << (a) << endl;
#define Print(x) cerr<<#x<<": "<<(x)<<endl
#define PrintVec(v) cerr<<#v<<": ";REP(__i,(v).size())cerr<<((v)[__i])<<", ";cerr<<endl
#endif


#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <climits>
#include <cfloat>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <set>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <sstream>
#include <complex>
#include <stack>
#include <queue>
#include <numeric>
#include <list>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <random>
#include <bitset>
#include <thread>
#include <chrono>
#include <vector>
#include <unordered_set>
#include <array>

using namespace std;
using namespace std::chrono;

typedef long long ll;
typedef unsigned long long ull;
template<class T> bool INRANGE(T x, T a, T b) { return a <= x && x <= b; }
template<class T> void amin(T& a, T v) { if (a > v) a = v; }
template<class T> void amax(T& a, T v) { if (a < v) a = v; }
template<class T> T CLAMP(T x, T a, T b) { return min(max(x, a), b); }
#define NG (-1)
#define BIG (static_cast<int>(1e9))
#define SZ(a) ((int)a.size()) 
#define SQ(a) ((a)*(a))
#define PI (acos(-1))
#define REP(i, n) for (int i = 0; (i) < (int)(n); ++ (i))
#define REP3(i, m, n) for (int i = (m); (i) < (int)(n); ++ (i))
#define REP_R(i, n) for (int i = (int)(n) - 1; (i) >= 0; -- (i))
#define REP3R(i, m, n) for (int i = (int)(n) - 1; (i) >= (int)(m); -- (i))
#define ALL(v) (v).begin(),(v).end()
#define ALL_OF(v, pred) all_of(ALL(v), [&](const auto& x) { return pred; })
#define ANY_OF(v, pred) any_of(ALL(v), [&](const auto& x) { return pred; })
#define NONE_OF(v, pred) none_of(ALL(v), [&](const auto& x) { return pred; })
#define COUNT_IF(v, pred) count_if(ALL(v), [&](const auto& x) { return pred; })

inline double RAD2DEG(double rad) { return rad * 180.0 / PI; }

static const double EPS = 1e-9;
inline int ROUND(double x) { return (int)(x + 0.5); }
inline bool ISINT(double x) { return fabs(ROUND(x) - x) <= EPS; }
inline bool ISEQUAL(double x, double y) { return fabs(x - y) <= EPS * max(1.0, max(fabs(x), fabs(y))); }
inline double SQSUM(double x, double y) { return x * x + y * y; }
typedef ll SCORE;
vector <SCORE> sScores;

class RandXor
{
public:
	RandXor()
	{
		init();
	}

	void init()
	{
		x = 123456789;
		y = 362436069;
		z = 521288629;
		w = 88675123;
	}

	void init(ll seed)
	{
		x = 123456789;
		y = static_cast<unsigned int>(362436069 + seed * 674890201);
		z = static_cast<unsigned int>(521288629 - seed * seed * 164907679);
		w = static_cast<unsigned int>(88675123 * seed);
	}


	inline unsigned int rand()
	{
		unsigned int t;
		t = (x ^ (x << 11)); x = y; y = z; z = w; return(w = (w ^ (w >> 19)) ^ (t ^ (t >> 8)));
	}

	unsigned int next() {
		return rand();
	}

	int nextInt(int num) {
		return rand() % num;
	}

	// [a,b)
	int next(int a, int b) {
		return a + (rand() % (b - a));
	}

	// [-d,d]
	int nextDelta(int d) {
		return next(-d, d + 1);
	}

	double nextD(double a, double b)
	{
		return a + (b - a) * nextDouble();
	}

	double nextDouble() {
		return (rand() + 0.5) * (1.0 / 4294967296.0);
	}

	template <class T>
	void randomShuffle(vector <T>& a)
	{
		const int n = SZ(a);
		for (int i = n - 1; i > 0; --i) {
			swap(a[i], a[nextInt(i + 1)]);
		}
	}

private:
	unsigned int x;
	unsigned int y;
	unsigned int z;
	unsigned int w;
};

// 秒を返す
inline double getTime()
{
	return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count() * 1e-6;
}

struct Pos
{
	short x;
	short y;

	Pos() = default;

	Pos(const Pos& other) {
		this->y = other.y;
		this->x = other.x;
	}

	Pos(short y, short x) {
		this->y = y;
		this->x = x;
	}

	bool operator< (const Pos& other) const
	{
		if (this->y < other.y)
		{
			return true;
		}
		else if (this->y == other.y)
		{
			if (this->x < other.x)
			{
				return true;
			}
		}

		return false;
	}

	Pos& operator=(const Pos& other)
	{
		this->y = other.y;
		this->x = other.x;
		return *this;
	}

	const Pos operator + (const Pos& other) const
	{
		return Pos(this->y + other.y, this->x + other.x);
	}

	Pos& operator += (const Pos& other)
	{
		this->y += other.y;
		this->x += other.x;

		return *this;
	}

	const Pos operator - (const Pos& other) const
	{
		return Pos(this->y - other.y, this->x - other.x);
	}

	Pos& operator -= (const Pos& other)
	{
		this->y -= other.y;
		this->x -= other.x;

		return *this;
	}

	bool operator== (const Pos& other) const
	{
		return this->y == other.y && this->x == other.x;
	}

	bool operator!= (const Pos& other) const
	{
		return !(*this == other);
	}

	bool isZero() const
	{
		return x == 0 && y == 0;
	}
};

template<class T>
vector<T> make_vector(int n, T t) {
	return vector<T>(n, t);
}

template<class ...Ts>
auto make_vector(int n, Ts ... ts) {
	return vector<decltype(make_vector(ts...))>(n, make_vector(ts...));
}

constexpr int NUM_DIR = 4;
constexpr int dy[NUM_DIR] = { 0, -1, 0, 1 }; // L,U,R,D
constexpr int dx[NUM_DIR] = { -1, 0, 1, 0 };
constexpr char dc[NUM_DIR] = { 'L', 'U', 'R', 'D'};


constexpr int N = 30;
constexpr int S = N*N;

const static int ROTATE[8] = {1, 2, 3, 0, 5, 4, 7, 6};
const static int TO[8][4] = {
	{1, 0, -1, -1},
	{3, -1, -1, 0},
	{-1, -1, 3, 2},
	{-1, 2, 1, -1},
	{1, 0, 3, 2},
	{3, 2, 1, 0},
	{2, -1, 0, -1},
	{-1, 3, -1, 1},
};
const pair<int, int> DIJ[4] = { {0, -1}, {-1, 0}, {0, 1}, {1, 0} };

vector < tuple <int, int, int> > getFullScoreRoute(const vector <vector <int> >& tiles, int sy = N / 2, int sx = N / 2) {

	vector <tuple <int, int, int> > tmp;
	bool used[N][N][4] = {};
	int y = sy;
	int x = sx;

	{
		int i = y;
		int j = x;
		int d = 0;
		if (TO[tiles[i][j]][d] != -1 && !used[i][j][d]) {
			// (si, sj) のタイルに sd 方向のタイルから侵入した状態からスタートして環状線の長さを求める
			int i2 = i;
			int j2 = j;
			int d2 = d;
			int length = 0;
			while (!used[i2][j2][d2]) {
				if (TO[tiles[i2][j2]][d2] == -1) {
					break;
				}
				length += 1;
				used[i2][j2][d2] = true;
				tmp.emplace_back(make_tuple(i2, j2, d2));

				d2 = TO[tiles[i2][j2]][d2]; // 次のタイルの方向
				used[i2][j2][d2] = true;
				tmp.emplace_back(make_tuple(i2, j2, d2));

				i2 += DIJ[d2].first;
				j2 += DIJ[d2].second;
				if (i2 >= N || j2 >= N || i2 < 0 || j2 < 0) { // 線路が途切れている
					break;
				}
				d2 = (d2 + 2) % 4;
			}

			y = i2;
			x = j2;
		}
	}

	return tmp;
}

pair <int, pair <int, int> > getFullScoreYX(const vector <vector <int> >& tiles, int sy = N / 2, int sx = N/2) {
	bool used[N][N][4] = {};
	int y = sy;
	int x = sx;
	int score = 0;

	{
		int i = y;
		int j = x;
		int d = 0;
		if (TO[tiles[i][j]][d] != -1 && !used[i][j][d]) {
			// (si, sj) のタイルに sd 方向のタイルから侵入した状態からスタートして環状線の長さを求める
			int i2 = i;
			int j2 = j;
			int d2 = d;
			int length = 0;
			while (!used[i2][j2][d2]) {
				if (TO[tiles[i2][j2]][d2] == -1) {
					break;
				}
				length += 1;
				used[i2][j2][d2] = true;
				d2 = TO[tiles[i2][j2]][d2]; // 次のタイルの方向
				used[i2][j2][d2] = true;
				i2 += DIJ[d2].first;
				j2 += DIJ[d2].second;
				if (i2 >= N || j2 >= N || i2 < 0 || j2 < 0) { // 線路が途切れている
					break;
				}
				d2 = (d2 + 2) % 4;
			}

			y = i2;
			x = j2;
			score = length;
		}
	}

	return make_pair(score, make_pair(y, x));
}

int getTwoCycleFullScore(const vector <vector <int> >& tiles, int sy = N / 2, int sx = N / 2) {

	vector <pair <int, vector <tuple <int, int, int> >> > ls;
	bool used[N][N][4] = {};

	auto cycle = make_vector(N, N, 4, make_pair(0, 0));

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int d = 0; d < NUM_DIR; d++) {
				if (TO[tiles[i][j]][d] != -1 && !used[i][j][d]) {
					int i2 = i;
					int j2 = j;
					int d2 = d;
					int length = 0;
					vector <tuple <int, int, int> > tmp;
					while (!used[i2][j2][d2]) {
						if (TO[tiles[i2][j2]][d2] == -1) {
							break;
						}
						length += 1;
						used[i2][j2][d2] = true;
						tmp.emplace_back(make_tuple(i2, j2, d2));
						d2 = TO[tiles[i2][j2]][d2];
						used[i2][j2][d2] = true;
						tmp.emplace_back(make_tuple(i2, j2, d2));
						i2 += DIJ[d2].first;
						j2 += DIJ[d2].second;
						if (i2 >= N || j2 >= N || i2 < 0 || j2 < 0) {
							break;
						}
						d2 = (d2 + 2) % 4;
					}
					if (i == i2 && j == j2 && d == d2) {
						ls.push_back(make_pair(length, tmp));

						for (const auto& tpl : tmp) {
							int tmpi;
							int tmpj;
							int tmpd;
							tie(tmpi, tmpj, tmpd) = tpl;
							cycle[tmpi][tmpj][tmpd].first = length;
						}
					}
				}
			}
		}
	}

	int score = 0;

	if (ls.size() <= 1) {
		score = 0;
	}
	else
	{
		sort(ALL(ls));
		score = ls[ls.size() - 1].first *ls[ls.size() - 2].first;
	};

	return score;
}

int getOneCycleFullScore(const vector <vector <int> >& tiles, int sy = N / 2, int sx = N / 2) {

	vector <pair <int, vector <tuple <int, int, int> >> > ls;
	bool used[N][N][4] = {};

	auto cycle = make_vector(N, N, 4, make_pair(0, 0));

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int d = 0; d < NUM_DIR; d++) {
				if (TO[tiles[i][j]][d] != -1 && !used[i][j][d]) {
					int i2 = i;
					int j2 = j;
					int d2 = d;
					int length = 0;
					vector <tuple <int, int, int> > tmp;
					while (!used[i2][j2][d2]) {
						if (TO[tiles[i2][j2]][d2] == -1) {
							break;
						}
						length += 1;
						used[i2][j2][d2] = true;
						tmp.emplace_back(make_tuple(i2, j2, d2));
						d2 = TO[tiles[i2][j2]][d2];
						used[i2][j2][d2] = true;
						tmp.emplace_back(make_tuple(i2, j2, d2));
						i2 += DIJ[d2].first;
						j2 += DIJ[d2].second;
						if (i2 >= N || j2 >= N || i2 < 0 || j2 < 0) {
							break;
						}
						d2 = (d2 + 2) % 4;
					}
					if (i == i2 && j == j2 && d == d2) {
						ls.push_back(make_pair(length, tmp));

						for (const auto& tpl : tmp) {
							int tmpi;
							int tmpj;
							int tmpd;
							tie(tmpi, tmpj, tmpd) = tpl;
							cycle[tmpi][tmpj][tmpd].first = length;
						}
					}
				}
			}
		}
	}

	int score = 0;

	if (ls.size() <= 0) {
		score = 0;
	}
	else
	{
		sort(ALL(ls));
		score = ls[ls.size() - 1].first;
	};

	return score;
}

void realMain(int caseID)
{
#ifndef SUBMISSION
	char input_path[256];
	char output_path[256];
	sprintf(input_path, "../../data/in/%04d.txt", caseID);
	sprintf(output_path, "../../data/out/output%04d.txt", caseID);
	ifstream fin(input_path);
	ofstream fout(output_path);
#endif

	RandXor randXor;
	vector <vector <int> > tiles;
	for (int i = 0; i < N; ++i)
	{
		string s;
		INPUT(s);

		vector <int> vi;
		for (char c : s)
		{
			vi.push_back(c - '0');
		}
		tiles.emplace_back(vi);
	}

	vector <vector <int> > firstTiles(tiles);


	vector <vector <int> > realBestTiles(tiles);
	int realBestScore = NG;


	const int NUM_MUL = 3;
	int sy = N / 2 + 1;
	int sx = N / 2 + 1;
	for (int multi=0;multi< NUM_MUL;multi++)
	{
		auto startClock = getTime();
		sx++;

		// 最高スコア
		{
			vector < vector <int> > state = tiles;
			vector < vector <int> > bestState = tiles;

			auto scyx = getFullScoreYX(state);
			auto score = scyx.first;
			auto lastYX = scyx.second;
			int bestScore = INT_MIN;

			const static double START_TEMP = 0.5; // 開始時の温度
			const static double END_TEMP = 0.5; // 終了時の温度
			const static double END_TIME = 1.5 / NUM_MUL; // 終了時間（秒）


			double temp = START_TEMP;   // 現在の温度

			long long steps;    // 試行回数
			for (steps = 0; ; steps++)
			{
				if (steps % 100 == 0)
				{
					const double time = getTime() - startClock;   // 経過時間（秒）
					if (time >= END_TIME)
					{
						break;
					}

					const double progressRatio = time / END_TIME;   // 進捗。開始時が0.0、終了時が1.0
					temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;
				}

				{
					// 近傍1 : ランダムな点を1つ回転
	//				const int y = CLAMP(lastYX.first  + randXor.nextDelta(1),0, N-1);
	//				const int x = CLAMP(lastYX.second + randXor.nextDelta(1),0, N-1);
					const int y = randXor.nextInt(N);
					const int x = randXor.nextInt(N);
					int& curD = state[y][x];
					const int bakD = curD;

					int deltaScore = 0;

					// 変更前のd日目のコンテストのスコアを引き、変更して、変更後のコンテストのスコアtを足す。
					deltaScore -= score;


					curD = ROTATE[curD];
					auto nextScYX = getFullScoreYX(state, sy, sx);
					deltaScore += nextScYX.first;
					const double probability = exp(deltaScore / temp); // 焼きなまし法の遷移確率


					if (steps % 10000 == 0)
					{
						//					cerr << " y=" << y << " x=" << x << endl;

						//					cerr << " deltaScore=" << deltaScore << " probability=" << probability << endl;
					}

					if (probability > randXor.nextDouble())
					{
						// 変更を受け入れる。スコアを更新
						score += deltaScore;
						lastYX = nextScYX.second;

						if (score > bestScore)
						{
							bestScore = score;
							bestState = state;
							//						cerr << " bestScore=" << bestScore << endl;

						}

					}
					else
					{
						// 変更を受け入れないので、元に戻す
						curD = bakD;
					}
				}
			}

			// 出力
			tiles = bestState;

			auto route = getFullScoreRoute(tiles, sy, sx);

			// ここから１番長いのを作成
			{
				int bestOneCycleScore = NG;
				auto bestTilesOneCycleScore = tiles;

				for (const auto& tpl : route)
				{
					int tmpi;
					int tmpj;
					int tmpd;
					tie(tmpi, tmpj, tmpd) = tpl;

					if (tiles[tmpi][tmpj] == 4 || tiles[tmpi][tmpj] == 5)
					{
						auto tmpTiles = tiles;

						tmpTiles[tmpi][tmpj] = ROTATE[tmpTiles[tmpi][tmpj]];


						const int sc = getOneCycleFullScore(tmpTiles, sy, sx);

						if (sc > bestOneCycleScore)
						{
							bestOneCycleScore = sc;
							bestTilesOneCycleScore = tmpTiles;
						}
					}
				}

				tiles = bestTilesOneCycleScore;
			}

			// ここから最高得点になるように作成
			{
				int bestTwoCycleScore = NG;
				auto bestTilesTwoCycleScore = tiles;

				for (const auto& tpl : route)
				{
					int tmpi;
					int tmpj;
					int tmpd;
					tie(tmpi, tmpj, tmpd) = tpl;

					if (tiles[tmpi][tmpj] == 4 || tiles[tmpi][tmpj] == 5)
					{
						auto tmpTiles = tiles;

						tmpTiles[tmpi][tmpj] = ROTATE[tmpTiles[tmpi][tmpj]];


						const int sc = getTwoCycleFullScore(tmpTiles, sy, sx);

						if (sc > bestTwoCycleScore)
						{
							bestTwoCycleScore = sc;
							bestTilesTwoCycleScore = tmpTiles;

						}
					}
				}

				tiles = bestTilesTwoCycleScore;

				cerr << " bestTwoCycleScore=" << bestTwoCycleScore << endl;

				if (bestTwoCycleScore > realBestScore)
				{
					realBestScore = bestTwoCycleScore;
					realBestTiles = tiles;
				}

			}
		}
	}


	sScores[caseID] = realBestScore;
	string ans;
	for (int y = 0; y < N; ++y)
	{
		for (int x = 0; x < N; ++x)
		{
			int tile = firstTiles[y][x];
			for (int rn = 0; rn < NUM_DIR; rn++)
			{
				if (tile == realBestTiles[y][x])
				{
					ans += string(1, '0' + rn);
					break;
				}
				tile = ROTATE[tile];
			}
		}
	}

	OUTPUT(ans);
	Print(ans);
}

vector <double> gCoefs;

int main(int argc, char* argv[])
{
	sScores.resize(1);
	if (argc > 1)
	{
		gCoefs.clear();
		gCoefs.resize(argc - 1);
		for (int i = 0; i < SZ(gCoefs); ++i)
		{
			gCoefs[i] = atof(argv[i + 1]);

		}
	}

#ifdef SUBMISSION
	realMain(0);
#else // SUBMISSION

	const int NUM_THREADS = 1;
	const int NUM_REPEATS = 10;
	const int ALL_TEST = NUM_THREADS * NUM_REPEATS;
	sScores.resize(ALL_TEST);

	SCORE totalScore = 0;

	for (int r = 0; r < NUM_REPEATS; ++r)
	{
		vector <std::thread*> p_threads(NUM_THREADS);
		for (int i = 0; i < NUM_THREADS; ++i)
		{
			p_threads[i] = new std::thread(&realMain, r * NUM_THREADS + i);
		}

		for (int i = 0; i < NUM_THREADS; ++i)
		{
			p_threads[i]->join();
		}

		for (int i = 0; i < NUM_THREADS; ++i)
		{
			delete p_threads[i];
		}
	}

	for (int i = 0; i < ALL_TEST; ++i)
	{
		cerr << i << " " << sScores[i] << endl;
		totalScore += sScores[i];
	}

	totalScore /= ALL_TEST;
	cout << "Average " << totalScore << endl; // optuna用
#endif // SUBMISSION

	return 0;
}


