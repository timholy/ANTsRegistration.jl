using Downloads

# Define the URL of the ANTs binary zip file
url = "https://github.com/ANTsX/ANTs/releases/download/v2.5.3/ants-2.5.3-ubuntu-22.04-X64-gcc.zip"

# Define paths for saving the zip file and extracting it
artifact_dir = joinpath(@__DIR__, "bin")  # Destination for binary files
zip_path = joinpath(artifact_dir, "ants-2.5.3-ubuntu-22.04-X64-gcc.zip")  # Using the correct file name

# Create the destination directory if it doesn't exist
isdir(artifact_dir) || mkpath(artifact_dir)

#Download the zip file
println("Downloading ANTs binary...")
Downloads.download(url, zip_path)

# Extract the downloaded zip file
println("Extracting ANTs binary...")

# Check the operating system for extraction command
if Sys.isunix() || Sys.isapple()
    # Use `unzip` command on Unix-like systems (Linux/macOS)
    run(`unzip -o $zip_path -d $artifact_dir`)
    elseif Sys.iswindows()
    # On Windows, use PowerShell to extract the zip file
    run(`powershell -Command "Expand-Archive -Path $zip_path -DestinationPath $artifact_dir"`)
else
    error("Unsupported operating system for automatic extraction.")
end

# Optionally remove the zip file after extraction
rm(zip_path)

println("ANTs binary installation completed.")
