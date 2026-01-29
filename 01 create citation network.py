"""
Note, this file is for reference only.
It would be used if you have manually downloaded results from the Web of Science, and want to turn them into a network.

If you're trying to replicate the work in our paper, skip this step and citations.txt will be used.
"""
from common import *

# load environment variables
load_dotenv()

BASE_DIR = os.environ.get("BASE_DIR")
print('base dir:', BASE_DIR)
BASE_DIR = Path(BASE_DIR)

fns = list(BASE_DIR.glob("**/*.txt"))
print(len(fns), 'files found')

work_range = 2010, 2019 # citing population

# helper for normalizing cited work strings
def normalize_c(x):
    a, *other = x.split(', ')
    aparts = a.split()
    if len(aparts) > 1:
        l, f, *_ = aparts
        f = f[0] # just first initial is fine!
        a = f"{l} {f}"

    # look for (v[0-9]+ and p[0-9]+) OR 'doi', to detect if it's a paper.
    # for books, we are just gonna throw away year
        
    if not (re.search(r'\b(v|p)[0-9]+\b', x) or 'doi' in x):
        try:
            int(other[0])
            other = other[1:]
        except:
            pass
    
    return ', '.join([a] + other)

cdict = {}

# collect citations
skipped_no_date = 0
for fn in tqdm(fns):
  with open(fn) as f:
    for li,l in enumerate(f):
      if li == 0: head = l.split(); continue
      parts = l.split('\t')
      py, doi = parts[head.index('PY')], parts[head.index('DI')]
      if py == '': 
        if doi != '':
          # get the part of the doi that looks like a date
          py = doi.split('.')
          py = [x for x in py if x.isdigit() and int(x) > 1900 and int(x) < 2020]
          if len(py):
            py = int(py[0])

          else:
            py = ''

      else:
        py = int(py)

      if py == '':
        # try using 66
        if parts[66] != '':
          py = int(parts[66].split('-')[0])

      if py == '':
        skipped_no_date += 1
        continue
      
      if py < work_range[0]-1 or py > work_range[1]+1:
        break # the whole file is not useful

      if py < work_range[0] or py > work_range[1]:
        continue

      citations = parts[29].split(';')
      citations = [x.strip().lower() for x in citations if x.strip() != '']
      citations = [normalize_c(x) for x in citations]
      citations = [x for x in citations if x != '[anonymous], communication']
      paper_id = parts[61]

      cdict[paper_id] = citations

# make sure each paper has a unique set of citations
for x in cdict:
  cdict[x] = list(set(cdict[x]))

index = {}
def i(name):
  if name not in index: index[name] = len(index)
  return index[name]

with open('networks/citations.txt', 'w') as outf:
  for f, cits in cdict.items():
    for t in cits:
      outf.write(f'{i(f)} {i(t)}\n')