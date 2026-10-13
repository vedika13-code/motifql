"""Assemble the MotifQL Playground: docs/index.html (served by GitHub Pages).

Run web/build-wasm.sh first to produce web/motifql.wasm, then:
    python3 web/build.py
"""
import base64, json, os
here = os.path.dirname(os.path.abspath(__file__))
t = open(os.path.join(here, 'template.html')).read()
w = base64.b64encode(open(os.path.join(here, 'motifql.wasm'), 'rb').read()).decode()
ex = json.dumps(json.load(open(os.path.join(here, 'examples.json')))).replace('</', '<\\/')
page = t.replace('"__WASM__"', '"' + w + '"', 1).replace('__EXAMPLES__', ex, 1)
head = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
        '<meta name="description" content="Write a DNA, RNA or protein motif, search FASTA files and see predicted chance hits. Runs in the browser.">\n'
        '<style>body{margin:0}img{max-width:100%}</style>\n</head>\n<body>\n')
out = os.path.join(here, '..', 'docs', 'index.html')
open(out, 'w').write(head + page + '\n</body>\n</html>\n')
print('wrote', os.path.normpath(out), os.path.getsize(out), 'bytes')
