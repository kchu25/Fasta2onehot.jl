function make_train_test_splits(num_reads; train_valid_test_ratio = [0.8, 0.1, 0.1])
    random_permuted_indices = shuffle(1:num_reads)
    train_end = floor(Int, num_reads * train_valid_test_ratio[1])
    valid_end = train_end + floor(Int, num_reads * train_valid_test_ratio[2])
    train_indices = random_permuted_indices[1:train_end]
    valid_indices = random_permuted_indices[train_end+1:valid_end]
    test_indices = random_permuted_indices[valid_end+1:end]
    return train_indices, valid_indices, test_indices, random_permuted_indices
end

function obtain_training_and_test_set(
    fastapath::String;
    train_valid_test_ratio = [0.9, 0.1, 0.0],
    float_type = Float32,
)
    onehotarr, onehotarr_shuffled = fasta2dummy(fastapath; F=float_type)
    
    train_indices, test_indices, _, random_permuted_indices = 
        make_train_test_splits(size(onehotarr, 4);
        train_valid_test_ratio = train_valid_test_ratio)

    training_set = @view onehotarr[:, :, :, train_indices]
    test_set = @view onehotarr[:, :, :, test_indices]
    training_set_shuffled = @view onehotarr_shuffled[:, :, :, train_indices]
    test_set_shuffled = @view onehotarr_shuffled[:, :, :, test_indices]
    return training_set, test_set, training_set_shuffled, test_set_shuffled
end
