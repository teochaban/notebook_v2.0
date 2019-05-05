int pai[MAXN], ranc[MAXN];

void init(int n){
    for(int i = 1; i <= n; i++){
        pai[i] = i;
    }
    memset(ranc, 1, sizeof(ranc));
}

int find(int i){
    if(pai[i] == i) return i;
    return pai[i] = find(pai[i]);
}

void join(int a, int b){
    a = find(a);
    b = find(b);
    if(a == b) return;
    if(ranc[a] > ranc[b]) ranc[a] += ranc[b], pai[b] = a;
    else ranc[b] += ranc[a], pai[a] = b;
}