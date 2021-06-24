#=
This module contains functions that parse and structure data into a
dict for use with the NGL Molecule Viewer component.
One or multiple input data files in the PDB or .cif.gz format can be
entered to return a dict to input as the `data` param of the component.
=#

"""
  single_split(string, sep)

Helper function to split `string` into 2 sub string based on `sep`
# Example
```julia
  a,b = DashBioUtils.single_split("hello.app",".")
```
"""
function single_split(string, sep)
  parts = split(string, sep)
  if length(parts) > 2
    ct = count(c -> c == only(sep), collect(string))
    throw("expected $(sep) once, found $(ct) in $(string)")
  end
  return parts
end

"""
get_highlights(string, sep, atom_indicator)

Helper function to set highlights using `string`, `sep`, `atom_indicator`

# Example
a,b = DashBioUtils.get_highlights("Hello.app", ".", "a")
"""
function get_highlights(string, sep, atom_indicator)
  residues_list = []
  atoms_list = []

  str_, _str = single_split(string, sep)

  for e in split(_str, ",")
    if occursin(atom_indicator, e)
      push!(atoms_list, e)
    else
      push!(residues_list, e)
    end
  end
  return (str_, Dict("atoms" => join(atoms_list, ","), "residues" => join(residues_list, ",")))
end

"""
get_data(data_path, pdb_id, color; reset_view=false, lc=true)
"""
function get_data(data_path, pdb_id, color; reset_view=false, loc=true)
  chain = "ALL"
  aa_range = "ALL"
  highlight_dic = Dict("atoms" => "", "residues" => "")

  # Check if only one chain should be shown
  if occursin(".", pdb_id)
    pdb_id, chain = single_split(pdb_id, ".")
    highlights_sep = "@"
    atom_indicator = "a"

    # Check if only a specified amino acids range should be shown:
    if occursin(":", chain)
      chain, aa_range = single_split(chain, ":")

      # Check if atoms should be highlighted
      if occursin(highlights_sep,aa_range)
        aa_range, highlight_dic = get_highlights(
          aa_range, highlights_sep, atom_indicator
        )
      end
    else
      if occursin(highlights_sep, chain)
        chain, highlight_dic = get_highlights(
          chain, highlights_sep, atom_indicator
        )
      end
    end  
  end
  if loc
    fname = [f for f in glob(string(data_path, pdb_id, ".*"))][1]
    if occursin("gz", fname)
      ext = split(fname, ".")[end-1]
      f = GZip.gzopen(fname,"r")
      rf = GZip.readline(f)
      content  = decode(rf, "UTF-8")
    else
      ext = split(fname, ".")[end]
      f = GZip.gzopen(fname,"r")
      content = GZip.readline(f)
    end
  else
    fname =  string(single_split(pdb_id, ".")[1], ".pdb")
    ext = split(fname, ".")[end]
    req = HTTP.request("GET", string(data_path, fname))
    content= decode(req.body, "UTF-8")
  end
  return Dict(
    "filename" => split(fname, "/")[end],
    "ext" => ext,
    "selectedValue" => pdb_id,
    "chain" => chain,
    "aaRange" => aa_range,
    "chosen" => highlight_dic,
    "color" => color,
    "config" => Dict("type" => "text/plain", "input" => content),
    "resetView" => reset_view,
    "uploaded" => false,
  )
end