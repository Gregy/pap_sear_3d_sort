/* 
 * File:   main.cpp
 * Author: gregy
 *
 * Created on October 10, 2013, 12:26 PM
 */

#include <cstdlib>
#include <vector>
#include <unordered_set>
#include <queue>
#include <string>
#include <set>
#include <bitset>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string.h>
#include <mpi.h>
#include <math.h>
// #define NDEBUG
#include <cassert>
#include <unistd.h>

using namespace std;

struct HEdge 
{
  unsigned int first;
  unsigned int second;
  unsigned int label;

public:
	HEdge(unsigned int first, unsigned int second, unsigned int label) {
		this->first = first;
		this->second = second;
		this->label = label;
	}
};
inline bool operator<(const HEdge& a, const HEdge& b)
{
	if(a.first < b.first) {
		return true;
	}
	if(a.first > b.first) {
		return false;
	}
	
	if(a.second <= b.second) {
		return true;
	}
	
	return false;
}


class Graph {
    vector<vector<unsigned int> > nodes;

public:
    Graph(unsigned int size):nodes(size) {
    }

    void addNeighbor(unsigned int node, unsigned int node2) {
        this->nodes.at(node).push_back(node2);
        this->nodes.at(node2).push_back(node);
    }

    vector<unsigned int> getNeighbors(unsigned int node) const {
        return this->nodes.at(node);
    }
    bool isNeighbor(unsigned int node, unsigned int node2) const{
        return find(this->nodes.at(node).begin(), this->nodes.at(node).end(), node2) != this->nodes.at(node).end();
    }
    void printGraph() {
        char *line = new char[this->nodes.size()+1];
        line[this->nodes.size()] = 0;

        for(unsigned int i=0;i<this->nodes.size();i++) {
            memset(line, '0', this->nodes.size());
            for(unsigned int ii=0;ii<this->nodes[i].size();ii++) {
                line[this->nodes[i][ii]] = '1';
            }

            cout << line << endl;
        }
        delete[] line;
    }
    void printToDot(bool printEnd = true) {
        cout << "graph graphname {" << endl;
        for(unsigned int i=0;i<this->nodes.size();i++) {
            for(unsigned int ii=0;ii<this->nodes[i].size();ii++) {
                if(this->nodes[i][ii] > i) {
		    cout << i << " -- " << this->nodes[i][ii] << ";" << endl;
                }
            }
        }
        if(printEnd) {
            cout << "}"<< endl;
        }
    }
    void printToDot(const vector<unsigned int> &path) {
        this->printToDot(false);
		double radius = path.size()*0.2;
        for(unsigned int i=0;i<path.size();i++) {
			double x = cos(((M_PI*2.0)/((double)path.size()))*((double)i))*radius;
			double y = sin(((M_PI*2.0)/((double)path.size()))*((double)i))*radius;
			cout << path[i] << "[label=\""<< path[i] << "\", pos=\"" << x << "," << y << "!\"];"<<endl;
			if(i != path.size()-1) {
				cout << path[i] << " -- " << path[i+1] << "[label=\"" << i << "\"];"<<endl;
			}
        }
    	cout << path.back() << " -- " << 0 << "[label=\"" << path.size()-1 << "\"];"<<endl;
        cout << "}" << endl;
    }

    unsigned int getNodeCount() const{
        return this->nodes.size();
    }
    
};

class FileGraphLoader {
public:
    static int load(string path, Graph** graph) {
        ifstream in_stream;
	string line;
	unsigned int size;
        
        try {
	in_stream.open(path);

        in_stream >> line;
	    size = std::stoi(line);
        }
        catch(exception& e) {
            return -1;
        }
	
        *graph = new Graph(size);
	for(unsigned int line_n=0;line_n<size; line_n++)
	{
            if(in_stream.eof()) {
                return -2;
            }
	    in_stream >> line;
            if(line.length() != size) {
                return -2;
            }
            for(unsigned int i =0;i<=line_n;i++) {
                if(line[i] == '1') {
                    (*graph)->addNeighbor(line_n, i);
                }
            }
	}

        return 0;
    }
};

class SearchState {
   vector<unsigned int> walkedPath;
   //state (next node to work on) for each node from walkedPath
   vector<unsigned int> nodeState;
   unordered_set<unsigned int> walkedNodes;

public:
    void pushNode(unsigned int node) {
        assert(this->walkedNodes.find(node)==this->walkedNodes.end());
        this->walkedPath.push_back(node);
        this->walkedNodes.insert(node);
        this->nodeState.push_back(0);
    }
    void removeLastNode() {
        this->walkedNodes.erase(this->walkedPath.back());
        this->nodeState.pop_back();
        this->walkedPath.pop_back();
    }
    unsigned int getState() {
        unsigned int ret = this->nodeState.back();
        return ret;
    }
    unsigned int getLastNode() {
        return this->walkedPath.back();
    }
    void setState(unsigned int s) {
        this->nodeState[this->nodeState.size()-1] = s;
    }
    bool isNodeInPath(unsigned int node) {
        if(this->walkedNodes.find(node) == this->walkedNodes.end()) {
            return false;
        }
		return true;
    }
    unsigned int getPathLength() {
        return this->walkedPath.size();
    }
    const vector<unsigned int> getPath() {
        return this->walkedPath;
    }
    
    void printPath() {
       const vector<unsigned int> path = this->getPath();
       for(unsigned int i=0;i<path.size();i++) {
           cout << path[i] << " ";
       }
       cout << endl;
    }

    unsigned int * asArray() {
        return this->walkedPath.data();
    }
};

class BFSStateGenerator {
    queue<SearchState*> states;
public:
    void generateStates(SearchState* initial, Graph * graph, unsigned int targetCount) {
		states.push(initial);
        while(states.size() < targetCount && states.size() > 0) {
            SearchState * cur = states.front();
            if(cur->getPathLength() == graph->getNodeCount()) {
                cout << "Prosel jsem cely graf v generujici fazi! Koncim s generovanim driv nez bys chtel!" << endl;
                break;
            }
            states.pop();
            const vector<unsigned int> neighbors = graph->getNeighbors(cur->getLastNode());
            for(unsigned int nextNode=cur->getState();nextNode<neighbors.size();nextNode++) {
                if(cur->isNodeInPath(neighbors[nextNode])) {
                    continue;
                }
                SearchState * next = new SearchState((*cur));
                next->pushNode(neighbors[nextNode]);
                states.push(next);
            }
            delete cur;
        }
    }
    void printStates() {
        for(unsigned int i=0; i<this->states.size();i++) {
            cout << i <<':';
			states.front()->printPath();
            states.push(states.front());
            states.pop();
		}
    }
    bool hasState() {
        return !this->states.empty();
    }
    SearchState * getNextState() {
		SearchState * cur = this->states.front();
        this->states.pop();
        return cur;
    }
    void clear() {
        while(!this->states.empty()) {
			this->states.pop();
        }
    }
	unsigned int getStatesCount() {
		return this->states.size();
	}
};

class DFSStateSolver {
    SearchState * state;
    const Graph * graph;

public:
    DFSStateSolver(SearchState * state, const Graph * graph) {
        this->state = state;
        this->graph = graph;
    }

    bool solve(const bool &breaker) {
        unsigned long uppermost = state->getLastNode();
		while(breaker) {
            //this->printPath();
            if(state->getPathLength() == graph->getNodeCount() && graph->isNeighbor(state->getLastNode(), 0)) {
                //found it!
                return true;
            }
            const vector<unsigned int> neighbors = graph->getNeighbors(state->getLastNode());
            if(state->getState() >= neighbors.size()) {
                if(state->getLastNode() == uppermost) {
                    return false;
                }
                state->removeLastNode();
                continue;
            }
            for(unsigned int nextNode=this->state->getState();nextNode<neighbors.size();nextNode++) {
                state->setState(nextNode+1);
                if(state->isNodeInPath(neighbors[nextNode])) {
                    continue;
                }
                state->pushNode(neighbors[nextNode]);
                break;
            }
        }
        return false;
    }
    void printPath() {
        state->printPath();
    }
    const vector<unsigned int> getPath() {
        return state->getPath();
    }
	set<HEdge> getPathVertices() {
       vector<unsigned int> path = state->getPath();
	   set<HEdge> s;

       for(unsigned int i=0;i<path.size()-1;i++) {
		   
		   s.insert(HEdge(path[i],path[i+1],i));
           cout << path[i] << " ";
       }
       cout << endl;
	   return s;
	}
	
};

class CLIArgumentsParser {
public:
    static char* getCmdOption(char ** begin, char ** end, const std::string & option)
	{
		char ** itr = std::find(begin, end, option);
		if (itr != end && ++itr != end)
		{
			return *itr;
		}
		return 0;
	}

	static bool cmdOptionExists(char** begin, char** end, const std::string& option)
	{
		return std::find(begin, end, option) != end;
	}
};

int main(int argc, char** argv) {

    const int STATES_PER_WORKER = 20;
    const int STATES_PER_THREAD= 15;

    char * filename = CLIArgumentsParser::getCmdOption(argv, argv + argc, "-f");

    if (!filename) {
        cout << "Chybi argument -f se souborem grafu" << endl;
		return -1;
    }

    Graph * graf;
    int res = FileGraphLoader::load(filename, &graf);
    if(res != 0) {
        cout << "chyba nacitani grafu ze souboru" << endl;
        return -1;
    }

	bool graphviz = CLIArgumentsParser::cmdOptionExists(argv, argv+argc, "-g");
    bool solved = false;

	MPI::Init();
    int groupSize = MPI::COMM_WORLD.Get_size();
    int myID = MPI::COMM_WORLD.Get_rank();
	
	if(CLIArgumentsParser::cmdOptionExists(argv, argv+argc, "-G")) {
		if(myID == 0)
			graf->printToDot();
		return 0;
	}

    if(groupSize < 2) {
        cout << "Tento program neumi bezet pouze s jednim nodem..." << endl;
        return -1;
    }

    
   	if(myID == 0) {
		cout << "Vypocetnich nodu:" << groupSize << endl;
		SearchState * state = new SearchState();
		state->pushNode(0);

		BFSStateGenerator * gen = new BFSStateGenerator();
		gen->generateStates(state, graf, groupSize*STATES_PER_WORKER);
		
        MPI::Status s;
        std::vector<bool> activeNodes(groupSize,true);
        bool answered = false;
        bool running= true;
        while(running) {
            MPI::COMM_WORLD.Probe(MPI::ANY_SOURCE, MPI::ANY_TAG, s);
            switch(s.Get_tag()) {
                case 10: //requesting work
                {
                    cout << "Node " << s.Get_source() << " pozaduje praci" << endl;
                    MPI::COMM_WORLD.Recv(NULL, 0,MPI::INT, s.Get_source(), s.Get_tag());
                    if(gen->hasState()) {
                        SearchState * stateToPass = gen->getNextState();
                        MPI::COMM_WORLD.Send(stateToPass->asArray(), stateToPass->getPathLength(), MPI::UNSIGNED, s.Get_source(), 10);
                        delete stateToPass;
                    }
                    else {
                        MPI::COMM_WORLD.Send(NULL, 0, MPI::UNSIGNED, s.Get_source(), 10);
                        activeNodes[s.Get_source()] = false;
                        activeNodes[0] = false;
                    }
                }
                    break;
                case 20: //ending computation
                {	
                    int count =s.Get_count(MPI::UNSIGNED);
                    if(count > 1 && !answered) { //first to arrive with data
                        answered = true;
                        activeNodes[0] = false;
                        activeNodes[s.Get_source()] = false;
						cout << "Node " << s.Get_source() << " nasel vysledek!" << endl;
						unsigned int * arr = new unsigned int[count];
						MPI::COMM_WORLD.Recv(arr, count, MPI::UNSIGNED, s.Get_source(), s.Get_tag());

						SearchState * state = new SearchState();
						for(int i=0;i<count;i++) {
							state->pushNode(arr[i]);
						}
                        delete[] arr;
						state->printPath();
                        solved = true;
                        if(graphviz) {
							graf->printToDot(state->getPath());
                        }
						cout << "Ukoncuji vypocet" << endl;
						gen->clear();
                        for(int i=1;i<groupSize;i++) {
                           if(!activeNodes[i]) {
                               continue;
                           }
                           MPI::COMM_WORLD.Send(NULL, 0, MPI::UNSIGNED, i, 20); 
                        }
                    }
                    else { //drop all others
						unsigned int * arr = new unsigned int[count];
                        MPI::COMM_WORLD.Recv(arr, count, MPI::UNSIGNED, s.Get_source(), s.Get_tag());
                        delete[] arr;
                        activeNodes[s.Get_source()] = false;
                    }
 				}                   

                    break;
                default:
                    cout << "Neznama zprava! Panika!" << endl;
            }
			
            running = false;
            for(unsigned int i=0;i<activeNodes.size();i++){
               if(activeNodes[i]) {
                   running = true;
               } 
            }
        }

		if(!solved && graphviz) {
			graf->printToDot();
		}
    } 
    else {//I am worker
        MPI::Status s;
		MPI::COMM_WORLD.Send(NULL, 0, MPI::UNSIGNED, 0, 10); //request work
		
        bool workerRun = true;
        while(workerRun) {
			MPI::COMM_WORLD.Probe(MPI::ANY_SOURCE, MPI::ANY_TAG, s);
			switch(s.Get_tag()) {
				case 10:
				{	
					int count =s.Get_count(MPI::UNSIGNED);
                    if(count < 1) {
                        workerRun = false;
                        break;
                    }
					unsigned int * arr = new unsigned int[count];
					MPI::COMM_WORLD.Recv(arr, count, MPI::UNSIGNED, s.Get_source(), s.Get_tag());

					SearchState * state = new SearchState();
					for(int i=0;i<count;i++) {
						state->pushNode(arr[i]);
					}
                    delete[] arr;
//                    {
//    int i = 7;
//    char hostname[256];
//    gethostname(hostname, sizeof(hostname));
//    printf("PID %d on %s ready for attach\n", getpid(), hostname);
//    fflush(stdout);
//    while (0 == i)
//        sleep(5);
//}

					cout << "Node " << myID << " prijal praci:";
					state->printPath();

                    BFSStateGenerator * gen = new BFSStateGenerator();
					gen->generateStates(state, graf, 4*STATES_PER_THREAD);
					cout << "Node " << myID << " nageneroval " << gen->getStatesCount() << " stavu pro lokalni vlakna." <<endl;
					bool askedForWork = false;

					#pragma omp parallel
                    {
                        while(workerRun) {
							SearchState * parState = NULL;
							#pragma omp critical(comm)
							{
								if(gen->hasState()) {
									parState = gen->getNextState();
								}
								else {
									if(!askedForWork) {
										askedForWork = true;
										MPI::COMM_WORLD.Send(NULL, 0, MPI::UNSIGNED, 0, 10); //request work
									}
								}
							}
							if(parState != NULL) {
								DFSStateSolver * solver = new DFSStateSolver(parState, graf);
								if(solver->solve(workerRun)) {
									#pragma omp critical(comm)
									{
										if(workerRun) {
											workerRun = false;
											MPI::COMM_WORLD.Send(solver->getPath().data(), solver->getPath().size(), MPI::UNSIGNED, 0, 20);
										}
									}
								}
								delete solver;
							}
                            else {
                                break;
                            }
                        }
                    }

                    delete gen;
				}
				break;
                case 20:
                {
					MPI::COMM_WORLD.Recv(NULL, 0, MPI::UNSIGNED, s.Get_source(), s.Get_tag());
					MPI::COMM_WORLD.Send(NULL, 0, MPI::UNSIGNED, 0, 20);
					workerRun = false;
                }
                break;
			}
        }
    }

	//cout << "Node " << myID << " konci." << endl;
    MPI::Finalize();


    return 0;
}