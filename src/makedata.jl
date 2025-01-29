function make_train_test_splits(num_reads; train_test_ratio = [0.8, 0.2])
    random_permuted_indices = shuffle(1:num_reads)
    train_end = floor(Int, num_reads * train_test_ratio[1])

    # map the indices of the shuffled data to the original data
    # shuffled index -> original index
    permute_map_train = @view random_permuted_indices[1:train_end]
    permute_map_test = @view random_permuted_indices[train_end+1:end]
    return permute_map_train, permute_map_test
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
    onehot_array::Array{T, 4}
    onehot_array_shuffled::Array{T, 4}
    training_set::AbstractArray{T, 4}
    training_set_shuffled::AbstractArray{T, 4}
    test_set::AbstractArray{T, 4}
    test_set_shuffled::AbstractArray{T, 4}
    permute_map_train::Vector{Int}
    permute_map_test::Vector{Int}
    data_seq_range::Vector{UnitRange{Int}}
    name::String
    cell_line::String
    source::String
    function onehot_data{T}(onehotarr, onehotarr_shuffled, 
            permute_map_train, permute_map_test, data_seq_range; 
            name="", cell_line="", source::String="") where {T <: Real}
        training_set          = @view onehotarr[:, :, :, permute_map_train]
        training_set_shuffled = @view onehotarr_shuffled[:, :, :, permute_map_train]
        test_set              = @view onehotarr[:, :, :, permute_map_test]
        test_set_shuffled     = @view onehotarr_shuffled[:, :, :, permute_map_test]

        new{T}(onehotarr, onehotarr_shuffled, 
               training_set, training_set_shuffled, 
               test_set, test_set_shuffled, 
               permute_map_train, permute_map_test, data_seq_range,
               name, cell_line,     source::String
               )
    end
end

get_num_seqs(data::onehot_data) = size(data.onehot_array, 4)
get_seq_len(data::onehot_data) = size(data.onehot_array, 2)

#=
Generate the meta data information (string)
for displaying the results in html
=#
function generate_meta_str(data::onehot_data)
    protein_name = isempty(data.name) ? "" : "<li> Protein name: $(data.name) </li>"
    cell_info = isempty(data.cell_line) ? "" : "<li> Cell line: $(data.cell_line) </li>"
    geo_accession = isempty(data.geo_accession) ? "" : "<li> GEO accession: $(data.geo_accession) </li>"
    num_seqs = "<li> Number of sequences: $(get_num_seqs(data)) </li>"
    seq_len = "<li> Sequence length: $(get_seq_len(data)) </li>"
    meta_str = protein_name * cell_info * geo_accession * num_seqs * seq_len
    return meta_str
end

# return how many rows we should have in the crosslink matrix
crosslink_mat_num_rows(data_instance::onehot_data) = 
    length(data_instance.data_seq_range)

function obtain_training_and_test_set(
    fastapaths::Vector{String};
    train_test_ratio = [0.9, 0.1],
    float_type = Float32,
    crosslink=false,
    name="",
    cell_line =""
)
    # load all the dna entries in the fasta files as onehot arrays
    onehotarr, onehotarr_shuffled, data_seq_range = 
        fasta2dummy(fastapaths; F=float_type, crosslink=crosslink)

    permute_map_train, permute_map_test = 
        make_train_test_splits(size(onehotarr, 4); train_test_ratio = train_test_ratio)
        
    return onehot_data{float_type}(onehotarr, 
                                   onehotarr_shuffled,
                                   permute_map_train, 
                                   permute_map_test,
                                   data_seq_range;
                                   name=name,
                                   cell_line=cell_line)
end

# for reproducing the same training and test set
struct data_instance_meta
    fastapaths::Vector{String}
    permute_map_train::Vector{Int}
    permute_map_test::Vector{Int}
end

fasta_get_meta_data(fasta_paths::Vector{String}, data_instance::onehot_data) = 
    data_instance_meta(fasta_paths, data_instance.permute_map_test, data_instance.permute_map_train)

# reproduce the same dataset using the saved meta data
function obtain_training_and_test_set(meta_data::data_instance_meta; float_type = Float32, name="", cell_line ="")
    onehotarr, onehotarr_shuffled, data_seq_range = 
        fasta2dummy(meta_data.fastapaths; F=float_type)
    permute_map_train, permute_map_test = 
        meta_data.permute_map_train, meta_data.permute_map_test
    return onehot_data{float_type}(onehotarr, 
                                   onehotarr_shuffled,
                                   permute_map_train, 
                                   permute_map_test,
                                   data_seq_range;
                                   name=name,
                                   cell_line=cell_line)
end


