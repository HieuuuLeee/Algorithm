
### 
```sh
long long C(int r, int n){
    if(n-r < r) r = n-r;
    if(r==0) p=1;
    while(r){
        p*=n; k*=r;
        int m = __gcd(p, k);
        p /= m; k/=m;
        n--; r--;
    }
    return p;
}
```