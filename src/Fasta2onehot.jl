module Fasta2onehot

using Random, SeqShuffle

include("dummy.jl")
include("makedata.jl")

export obtain_training_and_test_set, 
       which_dataset,
       crosslink_mat_num_rows

end
