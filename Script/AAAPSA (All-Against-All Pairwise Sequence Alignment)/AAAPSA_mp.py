# Input sequences
def needle_aln(seqleft, seqright, qinfo, sinfo, molecules):
    """
    cols = [ "OrthoGroup", "Qry.Acc", "Sbj.Acc", "Qry.Spec", "Sbj.Spec",
             "Pct.ID", "Pct.Sim", "Pct.Gaps", "nGaps", "Qry.Start",
             "Qry.End", "Sbj.Start", "Sbj.End", "Qry.Cov", "Sbj.Cov",
             "Score", "Program", "GapOpen", "GapExt.", "Matrix" ]
    """

    # All-against-all pairwise alignments ... <input_left> * <input_right> + 1 lines
    aln = psa.needle(moltype=molecules, qseq=seqleft, sseq=seqright)
    record = (f"{qinfo[0]}", f"{qinfo[2]}", f"{sinfo[2]}", f"{qinfo[1]}", f"{sinfo[1]}",
              f"{aln.pidentity}", f"{aln.psimilarity}", f"{aln.pgaps}", f"{aln.ngaps}", f"{aln.qstart}",
              f"{aln.qend}", f"{aln.sstart}", f"{aln.send}", f"{aln.query_coverage()}", f"{aln.subject_coverage()}",
              f"{aln.score:}", f"{aln.program}", f"{aln.gapopen}", f"{aln.gapextend}", f"{aln.matrix}")
    return "\t".join(record)


def water_aln(seqleft, seqright, qinfo, sinfo, molecules):
    """
    cols = [ "OrthoGroup", "Qry.Acc", "Sbj.Acc", "Qry.Spec", "Sbj.Spec",
             "Pct.ID", "Pct.Sim", "Pct.Gaps", "nGaps", "Qry.Start",
             "Qry.End", "Sbj.Start", "Sbj.End", "Qry.Cov", "Sbj.Cov",
             "Score", "Program", "GapOpen", "GapExt.", "Matrix" ]
    """

    # All-against-all pairwise alignments ... <input_left> * <input_right> + 1 lines
    aln = psa.needle(moltype=molecules, qseq=seqleft, sseq=seqright)
    record = (f"{qinfo[0]}", f"{qinfo[2]}", f"{sinfo[2]}", f"{qinfo[1]}", f"{sinfo[1]}",
              f"{aln.pidentity}", f"{aln.psimilarity}", f"{aln.pgaps}", f"{aln.ngaps}", f"{aln.qstart}",
              f"{aln.qend}", f"{aln.sstart}", f"{aln.send}", f"{aln.query_coverage()}", f"{aln.subject_coverage()}",
              f"{aln.score:}", f"{aln.program}", f"{aln.gapopen}", f"{aln.gapextend}", f"{aln.matrix}")
    return "\t".join(record)


if __name__ == "__main__":
    import sys
    import itertools
    import psa
    from multiprocessing import Pool
    from Bio import SeqIO

    if "p" in sys.argv[4].lower(): mol = "prot"
    else: mol = "nucl"

    itemsl = {seq.id: str(seq.seq) for seq in SeqIO.parse(sys.argv[1], "fasta")}
    if sys.argv[3] == "true": pass
    else: itemsr = {seq.id: str(seq.seq) for seq in SeqIO.parse(sys.argv[2], "fasta")}
    seql = []
    seqr = []
    qinfo = []
    sinfo = []
    molecules = []

    # sys.argv[1:2] = seqeucnes_left, right / sys.argv[3] = alignment_method / sys.argv[4] = molecular_type
    if sys.argv[3] == "true":
        for qid, sid in list(itertools.combinations(itemsl, r=2)):
            seql.append(itemsl[qid])
            seqr.append(itemsl[sid])
            qinfo.append(qid.split(":"))
            sinfo.append(sid.split(":"))
            molecules.append(mol)
    else:
        for qid, sid in list(itertools.product(itemsl, itemsr)):
            seql.append(itemsl[qid])
            seqr.append(itemsr[sid])
            qinfo.append(qid.split(":"))
            sinfo.append(sid.split(":"))
            molecules.append(mol)

    with Pool(processes=8) as pool:
        results = pool.starmap(needle_aln, zip(seql, seqr, qinfo, sinfo, molecules))

    for string in results: print(string)