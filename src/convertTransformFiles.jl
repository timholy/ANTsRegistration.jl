function Base.show(io::IO, s::ITKTransform)
# Define how to print the `id` field
    print(io, "ITKTransform(mode = $(s.mode), \n")
    print(io, "\tversion = $(s.version), \n")
    print(io, "\ttag = $(s.tag), \n")
    print(io, "\ttransform = $(s.transform), \n")
# Condense the description field to the first 50 characters (for example)
    tuple_len = length(s.parameters)
    if tuple_len > 10  # Limit the display to the first 10 elements
        short_parameters= s.parameters[1:10]  # Display the first 10 elements
        print(io, "\tdata = ($(join(short_parameters, ", ")), ... (truncated, total length: $tuple_len)), \n")
    else
        print(io, "\tdata = $(s.parameters) \n")  # Print the entire tuple if it's short
    end
    print(io, "\tfixedparameters = $(s.fixedparameters))")
end

"""ConvertTransformFile in ANTs"""
function convertTransformFile(tformfile::AbstractString, outputtxtfile::AbstractString)
    cmd = `ConvertTransformFile 2 $tformfile $outputtxtfile`
    run(cmd)
end

""" Load itk transform text file as a struct"""
function load_itktform(mode::AbstractString, tformtxtfile::AbstractString)
    if !isfile(tformtxtfile)
        @error "File does not exist"
    end
    lines = Vector{AbstractString}()
    open(tformtxtfile) do file
        for ln in eachline(file)
            push!(lines, ln)
        end
    end
    itktform_textlines2struct(mode, lines)
end
function load_itktform(tformtxtfile::AbstractString)
    if occursin("0GenericAffine.", tformtxtfile)
        mode = "GenericAffine"
    elseif occursin("1Warp.", tformtxtfile)
        mode = "Warp"
    elseif occursin("1InverseWarp.", tformtxtfile)
        mode = "InverseWarp"
    else
        @error "`mode` is unclear"
    end
    load_itktform(mode, tformtxtfile)
end

function itktform_textlines2struct(mode, lines)
    if isequal(lines[1], "#Insight Transform File V1.0") #Check file version
        ITKTransform(mode,
                     lines[1],
                     lines[2], 
                     lines[3][(findfirst(':', lines[3])+2):end],
                     Tuple(parse.(Float64, split(lines[4][(findfirst(':', lines[4])+2):end]))),
                     Tuple(parse.(Float64, split(lines[5][(findfirst(':', lines[5])+2):end]))))
    else
        @error "Transform file version does not match"
    end
end

function itktform_struct2textlines(tform::ITKTransform)
    lines = Vector{AbstractString}(undef, 5)
    lines[1] = tform.version
    lines[2] = tform.tag
    lines[3] = string("Transform: ", tform.transform)
    lines[4] = string("Parameters: ", join(tform.parameters, " "))
    lines[5] = string("FixedParameters: ", join(tform.fixedparameters, " "), "\n")
    return lines
end

""" Save as itk transform text file"""
function save_itktform(filename, tform::ITKTransform)
    lines = itktform_struct2textlines(tform)
    tformtxt = join(lines, "\n")
    open(filename, "w") do file
        write(file, tformtxt)
    end
end
