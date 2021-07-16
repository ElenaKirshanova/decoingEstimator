McEliece1 = {'n': 3488, 'k': 2720, 'w': 64}
McEliece2 = {'n': 4608, 'k': 3360, 'w': 96}
McEliece3 = {'n': 6960, 'k': 5413, 'w': 119}
McEliece4 = {'n': 8192, 'k': 6528, 'w': 128}


McElieceAll = [McEliece1, McEliece2, McEliece3, McEliece4]

# MMT ISD algorithm
# https://www.iacr.org/archive/asiacrypt2011/70730106/70730106.pdf
def MMT(params):

	n = params['n']
	k = params['k']
	w = params['w']

	p_step = 2
	l1_step = 2
	p = 10
	runtime_min = float("infinity")
	cost = {}

	while (p < round(w/3)):
		#log of representations
		# middle lists are m
		l2 = floor(log(binomial(p, round(p/2)),2))

		l1 = 2
		while(l1 < n - k - w + p - l2):
			l = l1+l2

			#list sizes
			S2 = binomial(round((k+l)/2), round(p/4)) #top-most
			S1 = S2^2 / 2^(l2) # middle
			S0 = S1^2 / 2^(l1) # bottom

			#probability of the correct weight distribution
			Prob = binomial((k+l)/2, p/2)^2 * binomial(n-k-l,w-p) / binomial(n,w)

			#runtime
			RTMMT = max(S2 , S1 , S0)/ Prob
			if(RTMMT < runtime_min):
				runtime_min = RTMMT
				T = [log(S2,2).n(), log(S1, 2).n(), log(S0, 2).n()]
				min_prob = log(Prob,2).n()
				mem = [log(S2,2).n(), log(S1,2).n()]
				p_min = p
				l_min = l
				minr = l2
				#print(p,l,runtime_min.n(),mem_min)
			l1 += l1_step
		p += p_step

	cost["runtime"] = log(runtime_min,2).n()
	cost["T"] = T
	cost["P"] = min_prob
	cost["mem"] = max(mem)
	cost["p"] = p_min
	cost["l"] = l_min
	return cost

#print('----------------MMT----------------')
#for param in McElieceAll:
#	res = MMT(param)
#	print('n =  ', param['n'], 'k = ', param['k'], 'w = ',  param['w'], res)



# BJMM ISD algorithm with depth 2 search tree
# https://www.iacr.org/archive/asiacrypt2011/70730106/70730106.pdf
def BJMM_d2(params):

	n = params['n']
	k = params['k']
	w = params['w']

	we = [0]*2
	runtime_min = float("infinity")
	cost = {}
	for eps in range(0, 12, 2): #number of additional 1's to represent 0+0

		p_step = 2
		l_step = 2
		p = 10

		while (p < round(w/3)):
			we[0] = p
			we[1] = round(we[0]/2)+eps

			l = 10
			while(l < n - k - w + p):

				#number of representations
				Rep = binomial(we[0],round(we[0]/2))*binomial(k+l-we[0], eps)
				r = floor(log(Rep, 2).n())


				#list sizes
				S2 = binomial(round((k+l)/2), round(we[1]/2))
				S1 = binomial(k+l, we[1]) / Rep # ==C2

				C2 = S2^2 / 2^r #expected number of matched pairs on lvl 1 (middle lvl)
				C1 = S1^2 / (2^(l - r)) #expected of matched pair on lvl 0 (bottom)

				#probability of the correct weight distribution
				Prob = binomial((k+l)/2, p/2)^2 * binomial(n-k-l,w-p) / binomial(n,w)
				RTBJMM = max(S2 , C2 , C1)/ Prob
				if(RTBJMM < runtime_min):
					runtime_min = RTBJMM
					T = [log(S2,2).n(), log(C2,2).n(), log(C1,2).n()]
					min_prob = log(Prob,2).n()
					Rep_min = [r]
					mem = [log(S2,2).n(), log(S1, 2).n()]
					p_min = p
					l_min = l
					eps_min = eps
					#print('new min:', log(runtime_min).n())
				l += l_step
			p += p_step

	cost["runtime"] = log(runtime_min,2).n()
	cost["T"] = T
	cost["P"] = min_prob
	cost["mem"] = max(mem)
	cost["p"] = p_min
	cost["l"] = l_min
	cost["eps"] = eps_min
	return cost


#print('----------------BJMM depth 2 ----------------')
#for param in McElieceAll:
#	res = BJMM_d2(param)
#	print('n =  ', param['n'], 'k = ', param['k'], 'w = ',  param['w'], res)

# BJMM ISD algorithm with depth 3 search tree
# https://www.iacr.org/archive/asiacrypt2011/70730106/70730106.pdf
def BJMM_d3(params):

	n = params['n']
	k = params['k']
	w = params['w']

	runtime_min = float("infinity")
	cost = {}
	p_step = 2
	l_step = 2
	we = [0]*3
	Rep = [0]*2
	r = [0]*2
	for eps1 in range(0, 12, 2): # additional 1's on lvl 1
		for eps2 in range(0, eps1,2): # additional 1's on lvl 2
			p = 10
			while (p < round(w/4)):

				we[0] = p
				we[1] = round(we[0]/2)+eps1
				we[2] = round(we[1]/2)+eps2

				l = 10
				while(l<n-k-w-p):

					#representtations on lvl 1
					Rep1 = binomial(we[0],round(we[0]/2))*binomial(k+l-we[0], eps1)
					r1 = floor(log(Rep1,2))

					#representtations on lvl 2
					Rep2 = binomial(we[1], round(we[1]/2))*binomial(k+l-we[1], eps2)
					r2 = floor(log(Rep2,2))

					S3 = binomial( round((k+l)/2), round(we[2]/2)) #list sizes on lvl 3
					S2 = binomial(k+l, we[2]) / Rep2 #list sizes on lvl 2
					S1 = binomial(k+l, we[1]) / Rep1 #list sizes on lvl 1


					C3 = S3^2 / (2^(r2)) #==S2
					C2 = S2^2 /(2^(r1-r2))
					C1 = S1^2 / (2^(l-r1))

					#probability of the correct weight distribution
					Prob = binomial((k+l)/2, p/2)^2 * binomial(n-k-l,w-p) / binomial(n,w)
					RTBJMM = max(S3 , C3 , C2, C1)/ Prob
					if(RTBJMM < runtime_min):
						runtime_min = RTBJMM
						mem_min = max(S3, S2, S1)
						p_min = p
						l_min = l
						eps_min = [eps1, eps2]
						#print('new min:', log(runtime_min).n())
					l += l_step
				p += p_step

			cost["runtime"] = log(runtime_min,2).n()
			cost["mem"] = log(mem_min).n()
			cost["p"] = p_min
			cost["l"] = l_min
			cost["eps"] = eps_min
	return cost

print('----------------BJMM depth 3 ----------------')
for param in McElieceAll:
	res = BJMM_d3(param)
	print('n =  ', param['n'], 'k = ', param['k'], 'w = ',  param['w'], res)
