function dna2dummy(dna_string::String, F::DataType)
    dummy = Dict(
        'A' => Array{F}([1, 0, 0, 0]),
        'C' => Array{F}([0, 1, 0, 0]),
        'G' => Array{F}([0, 0, 1, 0]),
        'T' => Array{F}([0, 0, 0, 1]),
        'N' => Array{F}([0, 0, 0, 0]),
        );
    v = Array{F,2}(undef, (4, length(dna_string)))
    found_n = false
    @inbounds for (index, alphabet) in enumerate(dna_string)
        alphabet_here = uppercase(alphabet)
        alphabet_here == 'N' && (found_n = true)
        v[:, index] = dummy[alphabet_here]
    end
    found_n && (@info "contains N in the string!")
    return v
end

# given an index find its order in the data_seq_range
function find_which(index, data_seq_range)
    for (ind, i) in enumerate(data_seq_range)
        if i.start ≤ index ≤ i.stop
            return ind
        end
    end
    error("index not found.")
end

function _dummy_(dna_reads, str_len, how_many_strings, F::DataType; data_seq_range=nothing)
    if isnothing(data_seq_range)
        _S_ = zeros(F, (4, str_len, 1, how_many_strings))
        @inbounds for i = 1:how_many_strings
            _S_[:, :, 1, i] = dna2dummy(dna_reads[i], F)
        end
        return _S_ # return a (4, str_len, 1, how_many_strings) array
    else
        crosslinked_position = str_len ÷ 2 + 1 # assume for now the center is the crosslinked position
        _S_ = zeros(F, (4+length(data_seq_range), str_len, 1, how_many_strings))
        @inbounds for i = 1:how_many_strings
            _S_[1:4, :, 1, i] = dna2dummy(dna_reads[i], F)
            which_crosslink = find_which(i, data_seq_range) + 4 # offset of 4 for the nucleotides
            _S_[which_crosslink, crosslinked_position, 1, i] = 1
        end
        return _S_ # return a (4+length(data_seq_range), str_len, 1, how_many_strings) array
    end
end

function data_2_dummy(dna_reads; F::DataType, data_seq_range=nothing)
    how_many_strings = length(dna_reads)
    @assert how_many_strings != 0 "There aren't DNA strings found in the input"
    _len_ = unique(length.(dna_reads))
    @assert length(_len_) == 1 "All DNA reads must be of same length"
    _S_ = _dummy_(dna_reads, _len_[1], how_many_strings, F; data_seq_range=data_seq_range)
    return _S_
end

#=
data_seq_range: Vector{UnitRange{Int}}
    length of data_seq_range is the number of datasets
    1st entry of data_seq_range is the range indices of the first dataset
    2nd entry of data_seq_range is the range indices of the second dataset
    ... and so on
=#
function get_fasta_range(seq_ends::Vector{Int})
    data_seq_range = Vector{UnitRange{Int}}(undef, length(seq_ends))
    data_seq_range[1] = 1:seq_ends[1]
    for i = 2:length(seq_ends)
        data_seq_range[i] = seq_ends[i-1]+1:seq_ends[i]
    end
    return data_seq_range
end

function which_dataset(data_seq_range::Vector{UnitRange{Int}}, index::Int)
    for i = 1:length(data_seq_range)
        index ∈ data_seq_range[i] && (return i)
    end
    return 0
end

function which_dataset(permute_map::Dict{Int, Int},
    data_seq_range::Vector{UnitRange{Int}}, index::Vector{Int})
    orig_ind = permute_map[index]
    return which_dataset(data_seq_range, orig_ind)
end

function fasta2dummy(fastapaths::Vector{String}; F::DataType, k = 1, crosslink = false)
    dna_heads = Vector{String}()
    dna_reads = Vector{String}()
    dna_reads_shuffled = Vector{String}()

    seq_counter = 0
    seq_ends = Int[]
    for fastapath in fastapaths
        f = open(fastapath)
        reads = read(f, String)
        close(f)
        for i in split(reads, ">")
            if !isempty(i)
                splits = split(i, "\n")
                this_read_head = splits[1]
                this_read = join(splits[2:end])
                push!(dna_heads, this_read_head)
                push!(dna_reads, this_read)
                push!(dna_reads_shuffled, seq_shuffle(this_read; k = k))
                seq_counter += 1
            end
        end
        push!(seq_ends, copy(seq_counter))
    end

    data_seq_range = get_fasta_range(seq_ends)

    # @info "data_seq_range: $data_seq_range"

    return data_2_dummy(dna_reads; F = F, 
                data_seq_range = crosslink ? data_seq_range : nothing), 
           data_2_dummy(dna_reads_shuffled; F = F, 
                data_seq_range = crosslink ? data_seq_range : nothing),
           data_seq_range
end
