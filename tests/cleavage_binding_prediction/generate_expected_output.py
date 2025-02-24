import pandas as pd
from crisprzip.kinetics import *

PRECISION = 1E-4

protospacer =  "AGACGCATAAAGATGAGACGCTGG"
targets = [protospacer,
           "AGTCGCATAAAGATGAGACGCGGG",
           "AGACCCATTAAGATGAGACGCGGG",
           "AGACGCATAACTATGAGACGCAGG",
           "AGACGCATAAAGATAAGCGGCCGG"]

protospacer2 =  "CAGTCATATCAGTCAGTACCTAGG"
targets2 = [protospacer2,
           "CAGTCATATCAGTCAGTACCTAGG",
           "CAGTCATATCAGTCAGTATCTTGG",
           "CAGTCATATCATTCACTACCTAGG",
           "CAATCATAACAATCAGTACCTCGG"]

protospacer3 =  "ACGTAGCTACTACATCGAGACTGG"
targets3 = [protospacer3,
           "ACGTAGCTACTACATCGAGACTGG",
           "ACTTAGGTACTACATCGAGACAGG",
           "ACGTAGCTACTACATCGATTCCGG",
           "ACGTAACTACCACATCGAGACTGG"]


time_pts = [10, 120, 3600]
k_bnd = [1E-2, 2E0, 5E2]

def average_params():

    avg_searcher = load_landscape("average_params")

    avg_df = pd.DataFrame(columns=['mm_pattern', 'exptype', 'k_on', 'time', 'fraction'])

    for tseq in targets:
        mmp = GuideTargetHybrid.from_cas9_offtarget(tseq, protospacer).get_mismatch_pattern()
        stc = avg_searcher.probe_target(mmp)
        for kb in k_bnd:
            for t in time_pts:
                f_clv = stc.get_cleaved_fraction(time=t, on_rate=kb).round(3)
                avg_df.loc[avg_df.shape[0]] = [str(mmp), 'cleavage', kb, t, f_clv]
        for kb in k_bnd:
            for t in time_pts:
                f_bnd = stc.get_bound_fraction(time=t, on_rate=kb).round(3)
                avg_df.loc[avg_df.shape[0]] = [str(mmp), 'binding', kb, t, f_bnd]

    return avg_df


def sequence_params():

    seq_searcher = load_landscape("sequence_params")

    seq_df = pd.DataFrame(columns=['protospacer', 'targetseq', 'exptype', 'k_on', 'time', 'fraction'])

    for tseq in targets:
        stc = seq_searcher.probe_sequence(protospacer, tseq)
        for kb in k_bnd:
            for t in time_pts:
                f_clv = stc.get_cleaved_fraction(time=t, on_rate=kb).round(3)
                seq_df.loc[seq_df.shape[0]] = [protospacer, tseq, 'cleavage', kb, t, f_clv]
        for kb in k_bnd:
            for t in time_pts:
                f_bnd = stc.get_bound_fraction(time=t, on_rate=kb).round(3)
                seq_df.loc[seq_df.shape[0]] = [protospacer, tseq, 'binding', kb, t, f_bnd]

    for tseq in targets2:
        stc = seq_searcher.probe_sequence(protospacer2, tseq)
        for kb in k_bnd:
            for t in time_pts:
                f_clv = stc.get_cleaved_fraction(time=t, on_rate=kb).round(3)
                seq_df.loc[seq_df.shape[0]] = [protospacer2, tseq, 'cleavage', kb, t, f_clv]
        for kb in k_bnd:
            for t in time_pts:
                f_bnd = stc.get_bound_fraction(time=t, on_rate=kb).round(3)
                seq_df.loc[seq_df.shape[0]] = [protospacer2, tseq, 'binding', kb, t, f_bnd]

    for tseq in targets3:
        stc = seq_searcher.probe_sequence(protospacer3, tseq)
        for kb in k_bnd:
            for t in time_pts:
                f_clv = stc.get_cleaved_fraction(time=t, on_rate=kb).round(3)
                seq_df.loc[seq_df.shape[0]] = [protospacer3, tseq, 'cleavage', kb, t, f_clv]
        for kb in k_bnd:
            for t in time_pts:
                f_bnd = stc.get_bound_fraction(time=t, on_rate=kb).round(3)
                seq_df.loc[seq_df.shape[0]] = [protospacer3, tseq, 'binding', kb, t, f_bnd]

    return seq_df


def main():
    avg_df = average_params()
    avg_df.to_csv("expected_output_avg.csv", index=False)

    seq_df = sequence_params()
    seq_df.to_csv("expected_output_seq.csv", index=False)


if __name__ == "__main__":
    main()
