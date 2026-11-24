import sys
S=[l.strip() for l in open(sys.argv[1])]
for l in sys.stdin:
    t=l.strip().split('_'); print('_'.join(t), ','.join(S[int(i)] for i in t), flush=False)
