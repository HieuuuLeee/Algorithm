# Find Bridge Offline on Graph
## Complexity: O(N + M)
## There are three case when you traversal the graph:
- visited[next_node] == false - the edge is part of DFS tree;
- visited[next_noe] == true && next_node ≠ previous_node - the edge is back edge to one of the ancestors;
- next_node == previous_node - the edge leads back to previous node in DFS tree.

## IMPLEMENT
### DECLARING VARIABLES
```sh
int n,m,timer,t1,t2;
vector<bool> visited;
vector<int> low,tin;
vector<vector<int>> next_v(111111);
vector<pair<int,int>> res;
```

### DFS FUNCTION O(N + M):
```sh
void dfs(int cur, int pre = -1){
	visited[cur] = true;
	tin[cur] = low[cur] = ++timer;

	for(int i : next_v[cur]){
		if(i == pre) continue;
		if(visited[i]) 
			low[cur] = min(low[cur], tin[i]);
		else{
			dfs(i, cur);
			low[cur] = min(low[cur], low[i]);
			if(low[i] > tin[cur])
				res.pb(mp(i,cur));
		}
	}
}
```

### FIND BRIDGE OFFLINE FROM ABITARY NODE ON GRAPH:
```sh
void find_bridges(){
	timer = 0;
	tin.assign(m,-1);
	low.assign(m,-1);
	visited.assign(m,false);

	for(int i=1; i<=m; i++){
		if(!visited[i]) 
			dfs(i);
	}
}
```
### MAIN FUNCTION
```sh
int main(){
	CURTIME();
    INFILE("in.txt");
    OUFILE("out.txt");

	cin>>n>>m; cout<<n<<"\n";
	For(i,0,n){
		cin>>t1>>t2;
		next_v[t1].pb(t2);
	}

	find_bridges();
	cout<<res.size();
	cout<<"\n";
	for(auto i : res)
		cout<<i.fi<<" "<<i.se<<"\n";
}
```