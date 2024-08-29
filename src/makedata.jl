function make_train_test_splits(num_reads; train_test_ratio = [0.8, 0.2])
    random_permuted_indices = shuffle(1:num_reads)
    train_end = floor(Int, num_reads * train_test_ratio[1])
    train_indices = random_permuted_indices[1:train_end]
    test_indices = random_permuted_indices[train_end+1:end]
    return train_indices, test_indices, random_permuted_indices
end

#=
permute_map: map the indices of the shuffled data to the original data, e.g.
 56  => 301, 35  => 134 says
    the 56th in the shuffled data is the 301st in the original data
    the 35th in the shuffled data is the 134th in the original data
data_seq_range: The range of the indices of the data in the original datasets, e.g.
     1:173
     174:399
     The first dataset consists of the sequences with indices 1 to 173
     The second dataset consists of the sequences with indices 174 to 399
=#
struct onehot_data{T}
    training_set::Array{T, 4}
    training_set_shuffled::Array{T, 4}
    test_set::Array{T, 4}
    test_set_shuffled::Array{T, 4}
    permute_map::Dict{Int, Int}
    data_seq_range::Vector{UnitRange{Int}}
    function onehot_data{T}(training_set, test_set, training_set_shuffled, test_set_shuffled, 
        permute_map, data_seq_range) where {T <: Real}
        new{T}(training_set, test_set, training_set_shuffled, test_set_shuffled, 
            permute_map, data_seq_range)
    end
end

# return how many rows we should have in the crosslink matrix
crosslink_mat_num_rows(data_instance::onehot_data) = 
    length(data_instance.data_seq_range)

function obtain_training_and_test_set(
    fastapaths::Vector{String};
    train_test_ratio = [0.9, 0.1],
    float_type = Float32,
)
    onehotarr, onehotarr_shuffled, data_seq_range = 
        fasta2dummy(fastapaths; F=float_type)

    train_indices, test_indices, random_permuted_indices = 
        make_train_test_splits(size(onehotarr, 4); train_test_ratio = train_test_ratio)

    # map the indices of the shuffled data to the original data
    # shuffled -> original
    permute_map = Dict(ind=>orig_ind for (ind, orig_ind) 
        in enumerate(random_permuted_indices))

    training_set = @view onehotarr[:, :, :, train_indices]
    test_set = @view onehotarr[:, :, :, test_indices]
    training_set_shuffled = @view onehotarr_shuffled[:, :, :, train_indices]
    test_set_shuffled = @view onehotarr_shuffled[:, :, :, test_indices]

    return onehot_data{float_type}(training_set, 
                                   training_set_shuffled,
                                   test_set, 
                                   test_set_shuffled, 
                                   permute_map, 
                                   data_seq_range)
end
