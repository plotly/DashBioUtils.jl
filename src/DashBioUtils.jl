module DashBioUtils

using Glob, GZip

using HTTP, StringEncodings, FASTX
using StatsBase

include("ngl_parser.jl")
include("protein_reader.jl")
include("xyz_reader.jl")

end # module
