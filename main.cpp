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
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string.h>
// #define NDEBUG
#include <cassert>

using namespace std;

class Graph {
    vector<vector<unsigned int>> nodes;

public:
    Graph(unsigned int size):nodes(size) {
    }

    void addNeighbor(unsigned int node, unsigned int node2) {
        this->nodes.at(node).push_back(node2);
        this->nodes.at(node2).push_back(node);
    }

    vector<unsigned int> getNeighbors(unsigned int node) {
        return this->nodes.at(node);
    }
    bool isNeighbor(unsigned int node, unsigned int node2) {
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
        for(unsigned int i=0;i<path.size()-1;i++) {
            cout << path[i] << " -- " << path[i+1] << "[label=\"" << i << "\"];"<<endl;
        }
    	cout << path.back() << " -- " << 0 << "[label=\"" << path.size()-1 << "\"];"<<endl;
        cout << "}" << endl;
    }

    unsigned int getNodeCount() {
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
        else {
            return true;
        }
    }
    unsigned int getPathLength() {
        return this->walkedPath.size();
    }
    const vector<unsigned int> getPath() {
        return this->walkedPath;
    }
};

class BFSStateGenerator {
    queue<SearchState*> states;
public:
    void generateStates(SearchState* initial, Graph * graph, unsigned int targetCount) {
	states.push(initial);
        while(states.size() < targetCount) {
            SearchState * cur = states.front();
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
};

class DFSStateSolver {
    SearchState * state;
    Graph * graph;

public:
    DFSStateSolver(SearchState * state, Graph * graph) {
        this->state = state;
        this->graph = graph;
    }

    bool solve() {
        unsigned long uppermost = state->getLastNode();
	while(true) {
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
    }
    void printPath() {
       const vector<unsigned int> path = state->getPath();
       for(unsigned int i=0;i<path.size();i++) {
           cout << path[i] << " ";
       }
       cout << endl;
    }
    const vector<unsigned int> getPath() {
        return state->getPath();
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

/*
 * 
 */
int main(int argc, char** argv) {

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

    SearchState * state = new SearchState();
    state->pushNode(0);

    DFSStateSolver * solver = new DFSStateSolver(state, graf);
    bool solved = solver->solve();
	bool graphviz = CLIArgumentsParser::cmdOptionExists(argv, argv+argc, "-g");
	if(graphviz) {
		cout << "/*"<< endl;
	}
    if(solved) {
        cout << "Hamiltonova cesta nalezena!" << endl;
        solver->printPath();
    }
    else {
		cout << "Hamiltonova cesta nenalezena :("<<endl;
    }
	if(graphviz) {
		cout << "*/" << endl;
		if(solved) {
			graf->printToDot(solver->getPath());
		}
		else{
			graf->printToDot();
		}
		cout.flush();
    }

    return 0;
}

