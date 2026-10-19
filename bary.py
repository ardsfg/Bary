import julia
julia.install()
jl = julia.Julia(compiled_modules=False)
from julia import Main

Main.include("/Users/ray/Documents/Geo/hello.jl")

# input = ""
# while(input != "quit"):
#     input = input()
    