module HDF5Helper
using HDF5

function create_file_group!(file::HDF5.File, group_name::String)
    create_group(file, group_name);
end

# creates 1D dataset
function create_dataset!(file::HDF5.File, group_name::String, dataset_name::String, datatype::DataType, chunk_size::Int)
    if group_name == ""
        create_dataset(file, dataset_name, datatype, ((0,),(-1,)), chunk = (chunk_size,));
    else
        create_dataset(file[group_name], dataset_name, datatype, ((0,),(-1,)), chunk = (chunk_size,));
    end
end

# appends 1D data to 1D dataset
function append_data!(file::HDF5.File, group_name::String, dataset_name::String, data::AbstractVector{Float64}, len_data::Int)
    if group_name == ""
        d = file[dataset_name];
    else
        d = file[group_name][dataset_name];
    end
    current_size = length(d);
    HDF5.set_extent_dims(d, (current_size+len_data,));
    d[current_size+1:current_size+len_data] = data;
end

# appends float to 1D dataset
function append_data!(file::HDF5.File, group_name::String, dataset_name::String, data::Float64)
    if group_name == ""
        d = file[dataset_name];
    else
        d = file[group_name][dataset_name];
    end
    current_size = length(d);
    HDF5.set_extent_dims(d, (current_size+1,));
    d[current_size+1] = data;
end

end

#= Example usage 

fname = tempname(); # temporary file

file = h5open(fname, "w")

group_name = "group1";
create_file_group!(file, group_name);

dataset_name = "dataset1";
chunk_size = 100;
create_dataset!(file, group_name, dataset_name, Float64, chunk_size);

chunk = rand(Float64, chunk_size)  # Example chunk data
@time append_data!(file, group_name, dataset_name, chunk, chunk_size)

file[group_name][dataset_name][101:end] == chunk

close(file)

=#