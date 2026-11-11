import collections, sys
def base(r):
    for m in ('_alignStrandInfo_','_AI_'):
        if m in r: r=r.split(m)[0]
    if '_split_' in r: r=r.split('_split_')[0]
    return r
d=sys.argv[1]
ind=collections.defaultdict(collections.Counter)
for line in open(d+'/results/read_classification.txt').readlines()[1:]:
    f=line.rstrip('\n').split('\t')
    ind[base(f[0]).split('_TRUTH_')[1]][f[4]] += 1
def hit(t):
    c=ind[t]; n=sum(c.values())
    good=c['Precise']+c['Precise_with_mutations']+c['Partial']+c['Partial_with_mutations']
    return n, good
for t in ['Replaced','Unedited','Del_only','Indel_8']:
    n,g=hit(t)
    if n: print(f'  {t:9s} n={n:5d}  called as carrying the edit: {g:5d} ({g*100/n:5.1f} %)')
