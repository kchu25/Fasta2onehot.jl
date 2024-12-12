module Fasta2onehot

using Random, SeqShuffle

include("dummy.jl")
include("makedata.jl")

export obtain_training_and_test_set, 
       onehot_data,
       which_dataset,
       crosslink_mat_num_rows,
       fasta_get_meta_data,
       get_num_seqs,
       get_seq_len

end
