hmmer_hits = set([])
hmmer_out_filenames = ["./tmp/fileKzsV9t"]

for hmmer_out_filename in hmmer_out_filenames:
    hmmer_out = open(hmmer_out_filename)
    for line in hmmer_out:
        line = line.strip()
        cols = line.split()
        print cols 
        if len(cols) < 14:
            continue
        
        sig = cols[13]
        if sig == "1":
            uscore = cols[0].rfind("_")
            if uscore != -1:
                hmmer_hits.add(cols[0][:uscore])
            else:
                print >> sys.stderr, "Warning: could not parse meta-assembly transcript id %s" % cols[0]
        
hmmer_hits_gene_ids = set([])
for t_id in hmmer_hits:
    dot = t_id.rfind(".")
    if dot != -1:
        g_id = t_id[:dot]
        hmmer_hits_gene_ids.add(g_id)
for t_id in hmmer_hits:
    dot = t_id.rfind(".")
    if dot != -1:
        g_id = t_id[:dot]
        if g_id in hmmer_hits_gene_ids:
            hmmer_hits.add(t_id)
for hit in hmmer_hits:
    print hit

