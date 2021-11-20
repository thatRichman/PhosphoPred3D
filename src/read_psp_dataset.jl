using DataFrames
using CSV

function read_psp_dataset(fh)
    psp_dat = CSV.File(fh;
        delim=",",
        header=1,
        threaded=true,
        tasks=4
    )
    psp_df = DataFrame(psp_dat)
    return psp_df
end

function read_pdb_dict()
    return CSV.File(PROJECT_DIR*"/data/PDB_IDS.tab", delim="\t") |> Dict
end