#include <bits/stdc++.h>
using namespace std;

#define ll long long

const ll MAX = 200001;
const ll MOD = 998244353;
const double PI = acosl(-1.0);

vector<ll> factorial(MAX), inverseFactorial(MAX);
ll limit = 1;
vector<ll> opt = { 0, 1 };


ll sub(ll x, ll y) {
    ll res = (x + MOD - y);
    if (res >= MOD) res -= MOD;
    return res;
}

ll power(ll x, ll y) {
    if (y == 0)
        return 1;
    if (y & 1)
        return (x * power(x, y - 1)) % MOD;
    ll res = power(x, y / 2);
    return (res * res) % MOD;
}

ll inverse(ll n) {
    return power(n, MOD - 2);
}

ll Binomial(ll n, ll r) {
    if (n < r)
        return 0;
    return (factorial[n] * ((inverseFactorial[r] * inverseFactorial[n - r]) % MOD)) % MOD;
}

void setup() {
    factorial[0] = factorial[1] = 1;
    for (int i = 1; i < MAX; i++)
        factorial[i] = (i * factorial[i - 1]) % MOD;

    inverseFactorial[MAX - 1] = inverse(factorial[MAX - 1]);
    for (ll i = MAX - 2; i >= 0; i--)
        inverseFactorial[i] = (inverseFactorial[i + 1] * (i + 1)) % MOD;
}

struct info {
    double x, y;
    info() {
        x = y = 0;
    }

    info(double x, double y): x(x), y(y) {}
};

inline info operator * (info a, info b) {
    return info(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}

inline info negation(info a) {
    return info(a.x, -a.y);
}

inline info operator-(info a, info b) {
    return info(a.x - b.x, a.y - b.y);
}

inline info operator+(info a, info b) {
    return info(a.x + b.x, a.y + b.y);
}

vector<info> alpha = {{ 0, 0 }, { 1, 0 }};

void checkInfo(ll num) {
    if (num <= limit) return;

    opt.resize(1 << num);

    for (ll i = 0; i < (1 << num); i++) {
        opt[i] = (opt[i >> 1] >> 1) + ((i & 1) << (num - 1));
    }

    alpha.resize(1 << num);

    while (limit < num) {
        double theta = 2 * PI / (1 << (limit + 1));
        for (ll i = 1 << (limit - 1); i < (1 << limit); i++) {
            alpha[i << 1] = alpha[i];
            double phi = theta * (2 * i + 1 - (1 << limit));
            alpha[(i << 1) + 1] = info(cos(phi), sin(phi));
        }

        limit++;
    }
}

void fastFourierTransform(vector<info> &list, ll n = -1) {
    if (n == -1) n = list.size();

    assert((n &(n - 1)) == 0);

    ll zc = __builtin_ctz(n);

    checkInfo(zc);

    ll rem = limit - zc;

    for (int i = 0; i < n; i++) {
        if (i < (opt[i] >> rem)) {
            swap(list[i], list[opt[i] >> rem]);
        }
    }

    for (int i = 1; i < n; i <<= 1) {
        for (int j = 0; j < n; j += 2 * i) {
            for (int k = 0; k < i; k++) {
                info ans = list[j + k + i] * alpha[k + i];
                list[j + k + i] = list[j + k] - ans;
                list[j + k] = list[j + k] + ans;
            }
        }
    }
}

vector<ll> product(vector<ll> &v1, vector<ll> &v2, ll equal = 0) {
    ll resSize = v1.size() + v2.size() - 1;
    int tmp = 0;
    while ((1 << tmp) < resSize) tmp++;

    checkInfo(tmp);

    vector<info> list1, list2;
    ll tpow = 1 << tmp;

    if (tpow > list1.size())
        list1.resize(tpow);

    for (int i = 0; i < v1.size(); i++) {
        ll val = (v1[i] % MOD + MOD) % MOD;
        list1[i] = info(val &((1 << 15) - 1), val>> 15);
    }

    fill(list1.begin() + v1.size(), list1.begin() + tpow, info{0, 0 });

    fastFourierTransform(list1, tpow);

    if (tpow > list2.size())
        list2.resize(tpow);

    if (equal)
        copy(list1.begin(), list1.begin() + tpow, list2.begin());
    else {
        for (int i = 0; i < v2.size(); i++) {
            ll val = (v2[i] % MOD + MOD) % MOD;
            list2[i] = info(val &((1 << 15) - 1), val>> 15);
        }

        fill(list2.begin() + v2.size(), list2.begin() + tpow, info{0, 0 });

        fastFourierTransform(list2, tpow);
    }

    double fraction = 0.25 / tpow;
    info root2(0, -1), root3(fraction, 0), root4(0, -fraction), root5(0, 1);

    for (int i = 0; i <= (tpow >> 1); i++) {
        int j = (tpow - i) &(tpow - 1);
        info info1 = (list1[i] + negation(list1[j]));
        info info2 = (list1[i] - negation(list1[j])) * root2;
        info info3 = (list2[i] + negation(list2[j])) * root3;
        info info4 = (list2[i] - negation(list2[j])) * root4;
        if (i != j) {
            info info5 = (list1[j] + negation(list1[i]));
            info info6 = (list1[j] - negation(list1[i])) * root2;
            info info7 = (list2[j] + negation(list2[i])) * root3;
            info info8 = (list2[j] - negation(list2[i])) * root4;
            list1[i] = info5 * info7 + info6 * info8 *  root5;
            list2[i] = info5 * info8 + info6 *  info7;
        }

        list1[j] = info1 * info3 + info2 * info4 *  root5;
        list2[j] = info1 * info4 + info2 *  info3;
    }

    fastFourierTransform(list1, tpow);
    fastFourierTransform(list2, tpow);

    vector<ll> finalRes(resSize);

    for (int i = 0; i < resSize; i++) {
        ll value1 = list1[i].x + 0.5;
        ll value2 = list2[i].x + 0.5;
        ll value3 = list1[i].y + 0.5;
        finalRes[i] = (value1 + ((value2 % MOD) << 15) + ((value3 % MOD) << 30)) % MOD;
    }

    return finalRes;
}

int main() {
    ios_base::sync_with_stdio(0);
    cin.tie(0);

    setup();

    int n;
    cin >> n;

    vector<vector < ll>> arr(n + 1, vector<ll> (32, 0));
    vector<ll> bits(32, 0);

    for (int i = 0; i < n; i++) {
        int x;
        cin >> x;
        int j = 0;
        while (x > 0) {
            if (x & 1) bits[j]++;
            x >>= 1;
            j++;
        }
    }

    for (int j = 0; j < 32; j++) {
        vector<ll> ones(bits[j] + 1, 0);
        vector<ll> zeroes(n - bits[j] + 1, 0);
        for (int i = 0; i <= n; i++) {
            if (i % 2 && i <= bits[j]) 
                ones[i] = Binomial(bits[j], i);
            if (i <= n - bits[j]) 
                zeroes[i] = Binomial(n - bits[j], i);
        }

        vector<ll> prod = product(ones, zeroes);
        for (int i = 1; i <= n; i++) {
            arr[i][j] = (arr[i - 1][j] + prod[i]) % MOD;
        }
    }

    int q;
    cin >> q;
    while (q--) {
        int qs;
        cin >> qs;

        ll ans = 0;

        for (int i = 0; i < 32; i++) {
            ans = (ans + ((arr[qs][i] * (1ll << i) % MOD) % MOD)) % MOD;
        }

        cout << ans << "\n";
    }

    return 0;
}
