#include <random>
#include <algorithm>
#include <iostream>
#include <string>
#include <fstream>
#include <chrono>
#include <vector>
#include <set>

#include "brkgaAPI/BRKGA.h"
#include "brkgaAPI/MTRand.h"
//#include "Node.h"

using namespace std;
using namespace std::chrono;

//typedef vector< Node > vNode;
typedef signed char**  sCC;

// Variávies Globais

double bestSolValue = 0;
vector<int> bestSolVet;
vector<set<int> > toolsPerTask;

// Dados do Problema

sCC Matrix_Graph; // Matriz de Chars - Mais rapido

int sizeOfMagazine = 1;
int numberOfTools = 1; // Linha
int numberOfTasks = 1; // Colunas

int FinalSolutionValue = INT_MAX;

unsigned K = 1;										// number of independent populations
unsigned MAXT = 4;									// number of threads for parallel decoding
unsigned X_INTVL = 100;								// exchange best individuals at every 100 generations
unsigned X_NUMBER = 2;								// exchange top 2 best
unsigned MAX_GENS = 100;							// run for 1000 gens
unsigned n = 10;										// size of chromosomes
unsigned p = 100;									// size of population
double pe = 0.20;									// fraction of population to be the elite-set
double pm = 0.10;									// fraction of population to be replaced by mutants
double rhoe = 0.70;									// probability that offspring inherit an allele from elite parent

vector <Subconjunto> Subconjuntos;

struct Subconjunto{  // struct para salvar as 2 posições para o swap
    int primeiraColuna;
    int segundaColuna;
};

// 1 - 2OPT // 2 - 2SWAP
const int localSeach = 1;

/*
 Quando rodar o IRACE utilizar uma seed FIXA
*/

//const long unsigned rngSeed = 2305843009213693951;	// seed to the random number generator

const long unsigned rngSeed = rand();


/* Funcoes */

double ktns(std::vector<int> taksOrder);
int getToolNeededSoonest(std::vector<int> spareTools,
                         std::vector<int> missingTasks);
std::vector<int> decodeChromosome(std::vector< double >& chromosome);


/*
 Reads the problem from a file specified by fileName
 */
void readProblem(string fileName)
{
  
  fileName = "Index/" + fileName;
  
  FILE* file = fopen (fileName.c_str(), "r");                         //input file
  
  if (!file) {
    //erro durante leitura
    exit(1);
  }
  
  fscanf(file,"%d",&numberOfTasks); // Num Ferrametnas - Coluna
  
  fscanf(file,"%d",&numberOfTools); // Num Ferrametnas - Linha
  
  fscanf(file,"%d",&sizeOfMagazine); // Tamanho Magazine
  
  Matrix_Graph = (signed char **) malloc (numberOfTools * sizeof (signed char*));
  
  for ( int k = 0; k < numberOfTools; k++){
    Matrix_Graph[k] = (signed char *) malloc (numberOfTasks * sizeof(signed char));
  }
  
  int num, i=0, j=0;
  
  while (!feof (file)){
    
    fscanf (file, "%d", &num);
    
    if(i == numberOfTools)
      break;
    
    Matrix_Graph[i][j] = num;
    
    j++;
    
    // reseta coluna
    if (j == numberOfTasks){
      i++;
      j=0;
    }
  }
  
  set<int> auxSet;
  for (int k = 0; k < numberOfTasks; k++){
    for(int j = 0; j < numberOfTools; j++){
      if((int)Matrix_Graph[j][k] == 1){
        auxSet.insert(j);
      }
    }
    toolsPerTask.push_back(auxSet);
    auxSet.clear();
  }
}

/*
 Generate Pairs of Combinations
*/

void Combinacao(int N, int K){
  string bitmask(K, 1); // K leading 1's
  bitmask.resize(N, 0); // N-K trailing 0's

Subconjunto sub;
vector<int> auxiliar;
// print integers and permute bitmask
do {
for (int i = 0; i < N; ++i) // [0..N-1] integers
{
if (bitmask[i]){
auxiliar.push_back(i);
}
}
sub.primeiraColuna = auxiliar[0];
sub.segundaColuna = auxiliar[1];
auxiliar.clear();
Subconjuntos.push_back(sub);
} while (std::prev_permutation(bitmask.begin(), bitmask.end()));
bitmask.clear();
}

void imprimeSubconjuntos(){
cout<< "Subconjuntos"<<endl;
for(int i = 0; i< Subconjuntos.size(); i++){
cout<< Subconjuntos[i].primeiraColuna<< ' '<<Subconjuntos[i].segundaColuna<<endl;
}
}

// Troca por intervalos
void dois_opt(){

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); // seed para evitar repetidos
vector<int> elementos;
int quantidade = Subconjuntos.size();

for (int i = 0; i <quantidade; i++){
elementos.push_back(i);
} 

shuffle (elementos.begin(), elementos.end(), std::default_random_engine(seed));
// SWAPLOCAL parametro
//int nTimes = round(quantidade * SWAP_LOCAL);

// fitness
int value = evaluationvalue[1];

for(int i = 0; i <10; i++){
vector<int>::iterator it = permutation.begin();
Subconjunto sub = Subconjuntos[elementos[i]];
reverse(it +(sub.primeiraColuna +1), it +(sub.segundaColuna+1));

// KTNS
consecutives();

// verifica ftiness
if (value <= evaluationvalue[1]){ // Piora de solucao, desfaz troca
reverse(it +(sub.primeiraColuna +1), it +(sub.segundaColuna+1));
evaluationvalue[1] = value;

} else {
value = evaluationvalue[1];                          
i = 0;
shuffle(elementos.begin(), elementos.end(), std::default_random_engine(seed));
}
}   
elementos.clear();
}

// Troca 2
void Pertubacao(){
  // seed para evitar repetidos
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); 
  vector<int> elementos;
  int quantidade = Subconjuntos.size();
  for (int i = 0; i <quantidade; i++){
    elementos.push_back(i);
  } 
  shuffle (elementos.begin(), elementos.end(), std::default_random_engine(seed));
   // quantas vezes serão executadas
  int nTimes = round(quantidade * SWAP_PERTUBATION);                        
  for(int i = 0; i < nTimes; i++){
    Subconjunto sub = Subconjuntos[elementos[i]];
    if (sub.primeiraColuna != -1 && sub.segundaColuna != -1){
      swap(permutation[sub.primeiraColuna], permutation[sub.segundaColuna]);
    }
  }
  // KTNS 
  consecutives(); 
  elementos.clear();  
}

/*
 Delta avaliacao
*/

int detaEvaluation(){
  return 0;
}

/*
 Busca Local  - Troca 2 opt
*/

void localSearch2Opt(int firstFitness,std::vector<int> chromossome){
  
  unsigned int index = 0;
  unsigned int size = int(chromossome.size());
  
  uint posA = 0;
  uint posB = 0;
  
  while (index < size) {
    
    
  }
  
}

/*
 Busca Local  - Troca 2 swap
*/

void localSearch2Swap(){

}

/*
 Refaz vetor solucao
*/


/*
 Metodo para refinar decodificacao
*/

int refineDecodeSolution(std::vector<int> firstPermutation, int firstFitness){
  
  
  // Escolhe qual busca local aplicar
  switch (localSeach) {
      
    case 1:
      localSearch2Opt(firstFitness,firstPermutation);
      break;
      
    case 2:
      localSearch2Opt(firstFitness,firstPermutation);
      break;
      
    default:
      localSearch2Opt(firstFitness,firstPermutation);
      break;
  }
  
  
  return 0;
}

/*
 Decoder Class from BRKGA
*/

std::vector<int> decodeChromosome(std::vector< double > chromosome){
  
  typedef std::pair< double, unsigned > ValueKeyPair;
  std::vector< ValueKeyPair > rank(chromosome.size());
  
  for(int i = 0; i < (int)chromosome.size(); i++) {
    rank[i] = ValueKeyPair(chromosome[i], i);
  }
  
  std::sort(rank.begin(), rank.end());
  
  std::vector< int > permutation;
  for(std::vector< ValueKeyPair >::const_iterator i = rank.begin(); i != rank.end(); ++i) {
    permutation.push_back(i->second);
  }
  
  return permutation;
}

class Decoder {
public:
  double decode(const std::vector< double >& chromosome) const  {
    
    double myFitness = 0.0;
    
    std::vector< int > permutation = decodeChromosome(chromosome);
    
    myFitness = ktns(permutation);
    
    // busca local
    
    myFitness = refineDecodeSolution(permutation,myFitness);
    
    return myFitness;
  }
};

/*
 KTNS function to evaluate fitness
 */

double ktns(std::vector<int> taksOrder){
  
  std::set<int> magazine;
  
  // Adiciona todas as ferramentas da tarefa 1 a caixa
  double fitness = toolsPerTask[taksOrder[0]].size();
  
  magazine = toolsPerTask[taksOrder[0]];
  
  for(int task = 1; task < numberOfTasks; task++){
    
    int emptySpaceInMag = sizeOfMagazine - (int)magazine.size();
    
    std::set<int> setTools = std::set<int>(toolsPerTask[taksOrder[task]].begin(),toolsPerTask[taksOrder[task]].end());
    
    //verifica se existe espaco vazio no magazine
    if(emptySpaceInMag >= setTools.size()){
      magazine.insert(setTools.begin(),setTools.end());
      fitness += setTools.size();
      
    }else{
      
      // E preciso fazer trocas
      
      std::vector<int>::iterator it;
      
      std::vector<int> diff(sizeOfMagazine,-1);
      std::vector<int> keepInMag(sizeOfMagazine,-1);
      std::vector<int> spareTools(sizeOfMagazine,-1);
      
      
      // The difference of two sets is formed by the elements that are present in the first set, but not in the second one.
      it = std::set_difference (setTools.begin(), setTools.end(),
                                magazine.begin(), magazine.end(),
                                diff.begin());
      
      diff.resize(it-diff.begin());
      
      it = std::set_difference (magazine.begin(), magazine.end(),
                                setTools.begin(), setTools.end(),
                                spareTools.begin());
      
      spareTools.resize(it-spareTools.begin());
      
      it = std::set_intersection (magazine.begin(), magazine.end(),
                                  setTools.begin(), setTools.end(),
                                  keepInMag.begin());
      
      keepInMag.resize(it-keepInMag.begin());
      
      
      if (diff.size() == 0){
        // diferenca entre o mag atual e as proxima tarefa for 0..pula tarefa
        break;
      }else{
        fitness += (int)diff.size();
        
        // diff contem as ferramentas presentes na tarefa corrente e nao estao no magazine
        // e preciso fazer diff.size trocas no magazine
        
        magazine = std::set<int>(keepInMag.begin(),
                                 keepInMag.end());
        
        if( (int)diff.size() == (sizeOfMagazine - (int)keepInMag.size()) ){
          // caso a diferenca seja igual ao tamanho do magazine..troca-se todas as pecas
          magazine.insert(diff.begin(),
                          diff.end());
          
        }else{
          // e preciso procurar qual ferramenta retirar do magazine
          for(int numTroca = 0; numTroca < (int)diff.size(); numTroca++){
            
            magazine.erase(getToolNeededSoonest(spareTools,
                                                std::vector<int>(taksOrder.begin() + task ,taksOrder.end())));
          }
          magazine.insert(diff.begin(),diff.end());
        }
      }
    }// Fim else
  }// Fim For
  
  return fitness;
}

int getToolNeededSoonest(std::vector<int> spareTools,
                         std::vector<int> missingTasks){
  
  typedef std::pair< int, int > Pair;
  std::vector< Pair > menosTrocas;
  
  for (int i = 0; i < (int)spareTools.size(); i++) {
    int tool = spareTools[i];
    int posNeeded = 0;
    while( posNeeded < missingTasks.size() && toolsPerTask[missingTasks[posNeeded]].find(tool) == toolsPerTask[missingTasks[posNeeded]].end()){
      posNeeded++;
    }
    menosTrocas.push_back(Pair(posNeeded,tool));
  }
  
  std::sort(menosTrocas.begin(),
            menosTrocas.end());
  
  
  return menosTrocas[0].second;
}



void printToolsVet(std::vector<int> vet){
  
  for (int i = 0; i < vet.size(); i++) {
    std::cout << vet[i] << " ";
  }
  std::cout<<std::endl;
  
}

// BRKGA

void setUpBRKGA(){
  
  const long unsigned rngSeed = 0;	// seed to the random number generator
  
  MTRand rng(rngSeed);				// initialize the random number generator
  
  Decoder decode;
  
  p = n * 100;
  
  // initialize the BRKGA-based heuristic
  BRKGA< Decoder, MTRand > algorithm(numberOfTasks,
                                     p,
                                     pe,
                                     pm,
                                     rhoe,
                                     decode,
                                     rng,
                                     K,
                                     MAXT);
  
  unsigned generation = 0;		// current generation
  
  do {
    algorithm.evolve();	// evolve the population for one generation
    
    if((++generation) % X_INTVL == 0) {
      algorithm.exchangeElite(X_NUMBER);	// exchange top individuals
    }
  } while (generation < MAX_GENS);
  
  bestSolValue = algorithm.getBestFitness();
  
  std::cout << "Best solution found has objective value = "
  << bestSolValue << std::endl;
  
  bestSolVet = decodeChromosome(algorithm.getBestChromosome());
  
}

/*
 Terminates all data structures.
 */
void termination()
{
  
  for (int i = 0;  i < numberOfTools; i++) {
    Matrix_Graph[i] = NULL;
  }
  
  Matrix_Graph = NULL;
  
  sizeOfMagazine = 0;
  numberOfTools = 0; // Linha
  numberOfTasks = 0; // Colunas
  bestSolVet.clear();
  bestSolValue = NULL;
  
  toolsPerTask.clear();
  
}

/*
 Prints the solution information to the output file
 */
void printSolution(string inputFileName, int solutionValue, double time, int run)
{
  string outputFileName =  "Results/Solution_Run_" + to_string(run) +"_" +  inputFileName ;
  
  ofstream fpSolution(outputFileName);				//file that contains the information about the solution of a problem instance
  
  fpSolution << "Instancia: " << inputFileName << std::endl;
  
  fpSolution << "Número de tarefas: " << numberOfTasks << std::endl;
  fpSolution << "Número de ferramentas: " << numberOfTools << std::endl;
  fpSolution << "Tamanho da caixa de ferramentas: " << sizeOfMagazine << std::endl;
  
  fpSolution << "Matriz de entrada:" << sizeOfMagazine << std::endl;
  for (int i = 0; i < numberOfTools; i++){
    for (int j = 0; j < numberOfTasks; j++){
      fpSolution << (int)Matrix_Graph[i][j] << " ";
    }
    fpSolution << endl;
  }
  
  fpSolution << "Run:" << run << std::endl;
  fpSolution << "Tempo decorrido:" << time << std::endl;
  fpSolution << "Melhor valor encontrado para as trocas:" << bestSolValue << std::endl;
  fpSolution << "Configuração da ordem das tarefas a serem executadas:"<< std::endl;
  
  for (int i = 0; i < numberOfTasks; i++){
    fpSolution << bestSolVet[i]+1 << " ";
  }
  fpSolution << endl;
  
  fpSolution << "Configuração da matriz final:"<< std::endl;
  
  for (int i = 0; i < numberOfTools; i++){
    for (int j = 0; j < numberOfTasks; j++){
      int index = bestSolVet[j];
      fpSolution << (int)Matrix_Graph[i][index] << " ";
    }
    fpSolution << endl;
  }
  
}


void multiRun(int *solutionValue, double *runningTime, string inputFileName, int run)
{
  
  //reads the problem data
  readProblem(inputFileName);
  
  //time taking
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  
  setUpBRKGA();
  
  //time taking
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  
  duration<double> time_span = duration_cast<duration<double> >(t2 - t1);
  
  //stores the solution value
  *solutionValue = bestSolValue;
  
  //stores the execution time
  *runningTime = time_span.count();
  
  //prints the solution to the file
  printSolution(inputFileName, *solutionValue, *runningTime, run);
  
  //terminates all data structures
  termination();
}
