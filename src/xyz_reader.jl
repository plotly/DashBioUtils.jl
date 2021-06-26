#=
XYZ reader
This module contains functions that can read an XYZ file and return a
Python dictionary with its contents.

=#

"""
read_xyz(datapath_or_datastring, is_datafile=true)

  Read data in .xyz format, from either a file or a raw string.
    `datapath_or_datastring`: Either the path to the XYZ file (can be relative
                                            or absolute), or a string corresponding to the content
                                            of an XYZ file (including newline characters).
    `is_datafile`: Either True (default) if passing the filepath to the data,
                                         or False if passing a string of raw data.
    `return`: A list of the atoms in the order that
                   they appear on the file, stored in
                   objects with keys "symbol", "x", "y",
                   and "z".
"""
function read_xyz(datapath_or_datastring; is_datafile=true)
  if is_datafile
    records = decode(read(datapath_or_datastring),"UTF-8")
  else
    req = Base.download(datapath_or_datastring)
    records = decode(read(req),"UTF-8")
  end 
  lines = split(records, "\n")

  atoms =Vector{Dict}(undef,0)

  for line in lines
    rs = match(
      r"^\s*([\w]+)\s+([\w\.\+\-]+)\s+([\w\.\+\-]+)\s+([\w\.\+\-]+)\s*", string(line)
    )
    if !(rs isa Nothing) && (length(rs.captures) == 4)
      atom = Dict(
        "symbol" => rs.captures[1],
        "x" => parse(Float64,rs.captures[2]),
        "y" => parse(Float64,rs.captures[3]),
        "z" => parse(Float64,rs.captures[4])
      )  
      push!(atoms, atom)  
    end
  end

  return atoms
end