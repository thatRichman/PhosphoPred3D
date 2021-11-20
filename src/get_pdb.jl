using CSV
using DataFrames
using HTTP
using DelimitedFiles
using BioStructures
using JSON
using Statistics
using Tar

include("../projvars.jl")
include("read_psp_dataset.jl")


function uniprot_acc_2_pdb_id(up_acc::AbstractArray, write_out::Bool=false)
    # convert Uniprot Accessions to PDB Codes
    ##
    data_out = Dict()

    url = "https://www.uniprot.org/uploadlists/?"
    headers = ["application/x-www-form-urlencoded"]

    # chunk in order to prevent overloading server, getting blacklisted
    chunked = Base.Iterators.partition(up_acc,400)
    for part in chunked
        params = Dict(
            "from" => "ACC+ID",
            "to" => "PDB_ID",
            "format" => "tab",
            "query" => join(part," ")
        )
        body = HTTP.escapeuri(params)
        resp = HTTP.request("GET", url*body, redirect=true)
        data = String(resp.body)
        data_splt = split(data, "\n")
        data_fmt = map(x -> split(x, "\t"), data_splt)[2:end-1] |> x->reduce(hcat,x)
        data_dict = Dict(zip(data_fmt[2,:], data_fmt[1,:]))
        merge!(data_out, data_dict)
        sleep(3)
    end

    if write_out
        open(PROJECT_DIR*"/data/PDB_IDS.tab", "w") do io
            write(io, "PDB_ID\tUP_ACC\n")
            for (k,v) in data_out
                write(io, k,"\t",v,"\n")
            end
        end
    end

    return data_out
end

function uniprot_acc_2_pdb_id(up_acc::AbstractString)
    params = Dict(
            "from" => "ACC+ID",
            "to" => "PDB_ID",
            "format" => "tab",
            "query" => up_acc,
        )
    url = "https://www.uniprot.org/uploadlists/?"
    body = HTTP.escapeuri(params)
    resp = HTTP.request("GET", url*body, redirect=true)
    data = String(resp.body)
    # data is returned as a tab-delimited string so we split on newlines 
    # then on tabs,
    # remove leading column header and trailing empty,
    # use reduce with hcat to convert the array of arrays to a matrix, 
    # return the second row
    data_splt = split(data, "\n")
    data_out = map(x -> split(x, "\t"), data_splt)[2:end-1] |> x->reduce(hcat,x)[2,:]
    return String.(data_out)
end
##

# fetch PDB structure files
# In the end we only want one structure per unique UP entry
# so we will select that with the highest quality metrics

"""
    get_quality_metrics_pdb(ids::AbstractArray)

    POST a GraphQL query to the RCSB PDB and return:
        title,
        absolute_percentile_clashscore,
        absolute_percentile_percent_ramachandran_outliers,
        absolute_percentile_percent_RSRZ_outliers,
        PDB_resolution,
"""
function get_quality_metrics_pdb(ids::AbstractArray)
    pdb_gql_endpoint = "https://data.rcsb.org/graphql?"
    join_ids = join(["\"$x\"" for x in ids], ",")
    query = """
        query getMetrics(\$pdb_ids: [String!]!) {
            entries(entry_ids: \$pdb_ids)
            {
                rcsb_id
                struct
                {
                    title
                }
                pdbx_vrpt_summary
                {
                    absolute_percentile_clashscore
                    absolute_percentile_percent_ramachandran_outliers
                    absolute_percentile_percent_RSRZ_outliers
                    PDB_resolution
                }
            }
        }
    """
    vars = """
    {
        "pdb_ids": [$join_ids]
    }       
    """
    # TODO once Diana.jl supports lists can transition
    params = Dict(
        "query" => query,
        "variables" => vars
    )
    resp = HTTP.request(
        "POST", 
        pdb_gql_endpoint, 
        ["Content-Type" => "application/json"], 
        JSON.json(params))
    return JSON.parse(String(resp.body))
end

"""
    find_highest_quality_struct(pdb_id, ap_ramachandran_outliers, app_rsrz_outliers, ap_clashscore, resolution)

    Given an Array of PDB IDs and arrays of associated quality metrics
    sort by resolution, then by median of absolute percentile quality scores
"""
function find_highest_quality_struct(pdb_id, ap_ramachandran_outliers, app_rsrz_outliers, ap_clashscore, resolution)
    
    vals = hcat(pdb_id, ap_ramachandran_outliers, app_rsrz_outliers, ap_clashscore, resolution)

    best_str = vals[1,:]
    for str in eachrow(vals[2:end,:])
        if str[5] < best_str[5] # resolution
            best_str = str
        elseif str[5] == best_str[5]
            if Statistics.median(str[2:4]) > Statistics.median(best_str[2:4])
                best_str = str
            end
        end
    end
    return (best_str[1])
end


"""
    sort_quality_metrics_up()

    Given a dict of pdb -> uniprot and quality metrics,
    subdivide quality metrics by associated uniprot ID
    and return a dataframe mapping uniprot IDs to the single "best" representative 3D structure
    by calling find_highest_quality_struct

"""
function sort_quality_metrics_up(map_dict::Dict, quality_metrics::Dict)
    
    rcsb_id = Vector{String}()
    up_id = Vector{String}()
    title = Vector{String}()
    appro =  Vector{Float64}()
    apprsrz = Vector{Float64}()
    apc = Vector{Float64}()
    res = Vector{Float64}()

    # create the DataFrame
    for entry in quality_metrics["data"]["entries"]
        push!(rcsb_id, entry["rcsb_id"])
        push!(up_id, map_dict[entry["rcsb_id"]])
        push!(title, entry["struct"]["title"])
        if entry["pdbx_vrpt_summary"] != nothing
            entry["pdbx_vrpt_summary"]["absolute_percentile_percent_ramachandran_outliers"] != nothing ? push!(appro, entry["pdbx_vrpt_summary"]["absolute_percentile_percent_ramachandran_outliers"]) : push!(appro, NaN)
            entry["pdbx_vrpt_summary"]["absolute_percentile_percent_RSRZ_outliers"] != nothing ? push!(apprsrz, entry["pdbx_vrpt_summary"]["absolute_percentile_percent_RSRZ_outliers"]) : push!(apprsrz, NaN)
            entry["pdbx_vrpt_summary"]["absolute_percentile_clashscore"] != nothing ? push!(apc, entry["pdbx_vrpt_summary"]["absolute_percentile_clashscore"]) : push!(apc, NaN)
            entry["pdbx_vrpt_summary"]["PDB_resolution"] != nothing ? push!(res, entry["pdbx_vrpt_summary"]["PDB_resolution"]) : push!(res, NaN)
        else
            push!(appro, NaN)
            push!(apprsrz, NaN)
            push!(apc, NaN)
            push!(res, NaN)
        end
    end

    quality_df = DataFrames.DataFrame(
        pdb_id = rcsb_id, 
        uniprot_id = up_id, 
        title = title, 
        ap_ramachandran_outliers=appro, 
        app_rsrz_outliers=apprsrz, 
        ap_clashscore=apc, 
        resolution=res)
    
    gdf = groupby(quality_df, "uniprot_id")
    best_struct_df = combine(gdf, [1,4,5,6,7] => ((p,a,b,c,d) -> (find_highest_quality_struct(p,a,b,c,d))) => :pdb_id)
    
    return best_struct_df
end


function main()   
    # convert uniprot IDs to PDB IDs
    if isfile(PROJECT_DIR*"/data/PDB_IDS.tab")
        @info "PDB_IDS.tab exists, reading..."
        pdb_dict = read_pdb_dict()
        psp_df = read_psp_dataset(PROJECT_DIR*"/data/Phosphorylation_site_dataset_psp.csv")
    else
        @info "parsing psp data and converting to PDB IDs"
        psp_df = read_psp_dataset(PROJECT_DIR*"/data/Phosphorylation_site_dataset_psp.csv")
        up_accs = unique(psp_df.ACC_ID)
        pdb_dict = uniprot_acc_2_pdb_id(up_accs, true)
    end

    # get quality metrics for all PDB IDs
    if isfile(PROJECT_DIR*"/data/quality_metrics.json")
        @info "quality metrics file found, reading..."
        f = open(PROJECT_DIR*"/data/quality_metrics.json", "r")
        quality_metrics_str = read(f, String)
        close(f)
        quality_metrics = JSON.parse(quality_metrics_str)
    else
        @info "Collecting quality metrics from PDB..."
        quality_metrics = pdb_dict |> values |> collect |> get_quality_metrics_pdb
        open(PROJECT_DIR*"/data/quality_metrics.json", "w") do f
            JSON.print(f, quality_metrics, 4)
        end
    end
    
    # select only the single best structure for each uniprot entry
    best_struct_df = sort_quality_metrics_up(pdb_dict, quality_metrics)

    joined_df = innerjoin(best_struct_df, psp_df, on = "uniprot_id" => "ACC_ID")
    CSV.write(PROJECT_DIR*"/data/best_structures_psp.csv", joined_df)

    # download the pdb files
    if !isfile(PROJECT_DIR*"/data/pdb_files/cif_files.tar.gz")
        downloadpdb(best_struct_df[!,2], dir=PROJECT_DIR*"/data/pdb_files/_cif", format=MMCIF)
        Tar.create(PROJECT_DIR*"/data/pdb_files/cif_files.tar.gz", pipeline(`gzip -9`, PROJECT_DIR*"/data/pdb_files/_cif"))
    else
        @warn "It looks like pdb files already exist, won't re-download"
    end
    return nothing
end


main()