## Big correlation
## Script by Julin Maloof
## Create Julia function to do correlation of large data sets in chunks.

## This uses an HDF5 disk-backed data array to keep memory usage reasonable.

## Currently the array will be overwritten everytime the function is called
## I need to deal with that at some point.

using HDF5 # Library for writing datasets to disk instead of keeping in memory

precision = Float32 # can use float 64 for things up to about 100000 columns

geno_array = rand(precision,200,200005) # create example data set

function chunkedcor(geno_array,chunk_size=NaN,verbose=true)
  # if no chunk size give the function will copmute on its own
  ncol = size(geno_array,2)

  if isnan(chunk_size)
    chunk_size = min(10000,floor(Int,ncol/10))
  end

  n_chunks = ceil(Int,ncol/chunk_size)

  if verbose
    println("chunk size: $chunk_size \n n.chunks: $n_chunks")
    println("creating results array")
  end

#for now results will get overwritten everytime the function is called.
  if exists(arrayfileid,"A")
    o_delete(arrayfileid,"A")
  end

  myarray = d_create(arrayfileid,
                      "A", #name of data object in file
                      datatype(precision),dataspace(ncol,ncol),
                        "chunk",(chunk_size,chunk_size),
                        "blosc",3) # blosc specifice compression. I am not sure about 3 but maybe that is compression level.

  if verbose
    println("array created")
  end

  for i in 1:n_chunks
    chunk_i_cols = ((i-1)*chunk_size + 1) : min(i*chunk_size,ncol) # so that we don't go over on the last column

    if verbose
      println("starting outer chunk i: $i")
    end

    for j in i:n_chunks
      chunk_j_cols = ((j-1)*chunk_size + 1) : min(j*chunk_size,ncol) # so that we don't go over on the last column
      tmp_cor = cor(geno_array[:,chunk_i_cols],geno_array[:,chunk_j_cols])
      if(i==j)
        tmp_cor = UpperTriangular(tmp_cor) # clean up if on diagonal
      end
      myarray[chunk_i_cols,chunk_j_cols] = tmp_cor
    end # for j
  end #for i
return(myarray)
end #chunkedcor

arrayfileid = h5open("cor_array.h5","w") #this is the file that will hold the correlation array

#test equality

nochunk_cor = UpperTriangular(cor(geno_array[:,1:1000])) #traditional correlation
nochunk_cor[1:5,1:5]

chunk_cor = chunkedcor(geno_array[:,1:1000])
chunk_cor[1:5,1:5]

chunk_cor = chunk_cor[:,:] #extract the array out of the object

isequal(chunk_cor,nochunk_cor) #no!

isequal(round(chunk_cor,2),round(nochunk_cor,2))
isequal(round(chunk_cor,3),round(nochunk_cor,3))

for i in 1:length(chunk_cor)
  if round(chunk_cor[i],3) != round(nochunk_cor[i],3)
    chunkval = chunk_cor[i]
    nochunkval = nochunk_cor[i] #because I can't figure out how to do variable sub with [] in a print command
    print("cell: $i chunked value: $chunkval nochunk value $nochunkval\n")
    #break
  end
end

# So things are not exactly the same but certainly good enough
# If Float64 is used then still not exactly the same but very, very close.

@time test1 = chunkedcor(geno_array[:,1:10003]) #1.6 seconds Float32; 3 sec Float64

@time test1 = chunkedcor(geno_array[:,1:100003]) # 228 sec Float32; 294 at float 64

@time test1 = chunkedcor(geno_array[:,1:200003]) # 656 sec Float32;

close(arrayfileid)
