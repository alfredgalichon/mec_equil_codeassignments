# Shortest Path Problem -- Bellman Operator
# Erica L Moszkowski

# Keys of graph are stop IDs
# Each entry of dict contains (destination, cost)
function load_graph(arc_file)
    graph = Dict{Int, Array{Tuple{Int, Float64}}}()

    infile = open(arc_file, "r")
    for line in readlines(infile)[2:end]
        elements    = split(line, ",")
        src_node    = parse(elements[1])
        dest_node   = parse(elements[2])
        weight      = parse(elements[3])

        entry       = (dest_node, float(weight))
        if src_node in keys(graph)
            push!(graph[src_node], entry)
        else
            graph[src_node] = [entry]
        end
    end
    close(infile)
    return graph
end

# # One step of Bellman operator
# function update_price(p, source, graph)
#     next_p  = Dict()
#     for node in keys(graph)
#         for neighbor
#         # if node == source
#         #     next_p[node] = 0.
#         # else
#             next_p[node] = min(p[node], minimum([cost + p[dest] for (dest, cost) in graph[node]]))
#         # end
#     end
#     return next_p
# end

# # Run Bellman til it converges
# function run_bellman(p, source, graph)
#     converged = false
#     # Run
#     i = 0
#     for i = 1:length(graph)
#         next_p = update_price(p, origin, graph)
#         if next_p == p
#             converged = true
#             break
#         else
#             p = next_p
#         end

#         if i == 10000
#             println("Maximum Number of Iterations Reached")
#         end
#     end

#     return p, i, converged
# end

function bellman_ford(graph, source)

    # Initialize
    distance    = Dict()
    predecessor = Dict()
    for node in keys(graph)
        distance[node]    = Inf
        predecessor[node] = -1
    end
    distance[source] = 0.

    # Relax edges
    for i in 1:length(graph)-1
        for node1 in keys(graph)
            for (node2, weight) in graph[node1]
                if distance[node1] + weight < distance[node2]
                    distance[node2]     = distance[node1] + weight
                    predecessor[node2]  = node1
                end
            end
        end
    end

    # Check for nonnegative weight cycle
    for node1 in keys(graph)
        for (node2, weight) in graph[node1]
            if distance[node1] + weight < distance[node2]
                info("Graph contains negative-weight cycle")
            end
        end

    end

    return distance, predecessor
end

function find_shortest_path(predecessors, source, dest)
    shortest_path = [dest]
    curr_node = dest
    while curr_node != source
        curr_node = predecessors[curr_node]
        push!(shortest_path, curr_node)
    end
    reverse(shortest_path)
end

# Initialize
graph   = load_graph("data/arcs.csv")
n_nodes = length(graph)

origin = 452
dest   = 471

distance, predecessor = bellman_ford(graph, origin)

# Find the shortest path
shortest_path, cost = find_shortest_path(predecessor, origin, dest)
