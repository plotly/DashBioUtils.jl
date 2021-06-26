

#=
This module includes functions that extract data from FASTA files into
dictionaries.
Attributes: _DATABASES (dict): A dictionary that translates the
    database sepecified in the description line of the FASTA file to
    its constituent metadata fields.
=#

# information on database header formats, taken from
# https://en.wikipedia.org/wiki/FASTA_format
_DATABASES = Dict(
    "gb" => ["accession", "locus"],
    "emb"=> ["accession", "locus"],
    "dbj" => ["accession", "locus"],
    "pir" => ["entry"],
    "prf" => ["name"],
    "sp" => ["accession", "entry name", "protein name", "organism name",
           "organism identifier", "gene name", "protein existence",
           "sequence version"],
    "tr" => ["accession", "entry name", "protein name", "organism name",
           "organism identifier", "gene name", "protein existence",
           "sequence version"],
    "pdb" => ["entry", "chain"],
    "pat" => ["country", "number"],
    "bbs" => ["number"],
    "gnl" => ["database", "identifier"],
    "ref" => ["accession", "locus"],
    "lcl" => ["identifier"],
    "nxp" => ["identifier", "gene name", "protein name", "isoform name"]
)

"""
  read_fasta(datapath_or_datastring; is_datafile=true)

Read a file in FASTA format, either from a file or from a string of raw
data.

`datapath_or_datastring`: Either the path to the FASTA file (can be relative
                                        or absolute), or a string corresponding to the content
                                        of a FASTA file (including newline characters).

`is_datafile`: Either `true` (default) if passing the filepath to the data,
                                     or `false` if passing a string of raw data.

`return`: A list of protein objects, each containing a
                     description (based on the header line) and the amino
                     acid sequence with, optionally, all non-amino-acid
                     letters removed.
"""
function read_fasta(datapath_or_datastring::String; is_datafile=true)
  raw_data = []
  # open file if given a path
  if is_datafile
    records = open(FASTA.Reader, datapath_or_datastring)
  else
    req = Base.download(datapath_or_datastring)
    records = open(FASTA.Reader, req)
  end 
  
  cl_records = collect(records)
  len_records =  length(collect(cl_records))
  fasta_data = Vector{Dict}(undef,len_records)
  for (idx, val) in enumerate(cl_records)
    dt = decode(val.data, "UTF-8")
   fasta_data[idx] = Dict("description" => decode_description(string(dt[val.identifier], dt[val.description])),
         "sequence" => string(dt[val.sequence]))
  end
  return fasta_data
end

"""
  decode_description(description)

Parse the first line of a FASTA file using the specifications of
    several different database headers (in _DATABASES).
    `description`: The header line with the initial `>``
                                 removed.
    `return`: A dictionary for which each key-value pair comprises
                   a property specified by the database used and the
                   value of that property given by the header. If the
                   database is not recognized, the keys are given as
                   'desc-n' where n is the position of the property.
"""
function decode_description(description)
  if description == ""
    return Dict("-1" => "no description")
  end
  decoded = Dict()
  desc = split(description,"|")
  if haskey(_DATABASES, desc[1])
    db_info = _DATABASES[desc[1]]
    if desc[1] in ["sp", "tr"]
      decoded["accession"] = desc[2]
      # using regex to get the other information      
      rs = match(
        r"([^\s]+)(.*)\ OS=(.*)\ OX=(.*)\ GN=(.*)\ PE=(.*)\ SV=(.*)$", desc[3]
      )
      for i in 3:length(db_info)
        decoded[db_info[i]] = rs.captures[i-2]
      end
    else
      # shift by one, since first section in header describes
      # the database
      for i in 1:length(desc)-1
        decoded[db_info[i]] = desc[i+1]
      end
    end
  else
    if length(desc) > 1
      for i in 1:length(desc)-1
        decoded[string(i)] = desc[i+1]
      end
    else
      ecoded["Header"] = desc[0]
    end
  end
  return decoded
end

