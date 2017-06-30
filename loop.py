N = 4
for i in range(N):
	for j in range(i+1,N,1):
		for u in range(i,N,1):
			for v in range(u+1, N, 1):
				if j < v or i < u:
					print (i+1,j+1),(u+1,v+1)