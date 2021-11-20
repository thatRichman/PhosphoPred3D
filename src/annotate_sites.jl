using BioStructures
using BioSequences
using BioAlignments
include("../projvars.jl")
include("read_psp_dataset.jl")

"""
    find_site(query::LongAminoAcidSeq, target::LongAminoAcidSeq, approx_max_dev::Int8=2)

    Attempt to find the query sequence in the target, first by exact match and then approximate
    
    approx_max_dev controls the maximum allowable sequence levenshtein distance, default 2.

"""
function find_site(query::Union{AminoAcidSeq,LongAminoAcidSeq}, target::Union{AminoAcidSeq,LongAminoAcidSeq}, approx_max_err::Int64=2)
    query_pp = ExactSearchQuery(query)
    hit = findfirst(query_pp, target)
    if !isnothing(hit)
        return (hit,0)
    else
        # attempt approximate search
        query_pp = ApproximateSearchQuery(query)
        approx_hit = approxsearch(target, query_pp, approx_max_err)
    end
    if approx_hit != 0:-1
        return (approx_hit, BioSequences.levenshtein_distance(query, target))
    else
        return nothing
    end
end

function preprocess_site_sequences(seqs::Array{String})
    return LongAminoAcidSeq.(replace.(seqs, "_" => ""))
end

function main()
    psp_df = DataFrame(CSV.File(PROJECT_DIR*"/data/best_structures_psp.csv"))
    psp_df[!,:pdb_id ] = String.(psp_df[!,:pdb_id])
    psp_site_pairs = psp_df[!,["SITE_+/-7_AA","pdb_id"]]
    psp_site_pairs[!,"SITE_+/-7_AA"] = preprocess_site_sequences(psp_site_pairs[!,"SITE_+/-7_AA"])

    sites = []
    for row in eachrow(psp_site_pairs)
        struc = read(PROJECT_DIR*"/data/pdb_files/cif_files/$(row[:pdb_id]).cif", MMCIF)
        for mod in struc
            for chain in mod
                found_site = find_site(LongAminoAcidSeq(row["SITE_+/-7_AA"]), LongAminoAcidSeq(chain, standardselector))
                if !isnothing(found_site)
                    push!(sites, tuple(row["pdb_id"], row["SITE_+/-7_AA"], mod, chain, found_site))
                end
            end
        end
    end
    print(sites)
end

main()