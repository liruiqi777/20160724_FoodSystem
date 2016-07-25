/*****************************************************************************
 File:   foodSystems.cc
 Author: Adela Yang, Ruiqi (Sally) Li, Lucy Skinner
         Bowdoin College '16/'18/'16
 Date:   May 2016
 
 Description:
 Implement a food systems demo using matrixes
 
 Running instructions:
 save folder on Desktop
 open Terminal application
 type in:   cd Desktop/foodSystems
 press enter
 then type in: g++ -Wall -o food foodSystems.cc
 press enter
 then type in: ./food
 press enter
 
 NOTES:
 Matrix A = refers to the big first matrix (influence)
 Matrix B = refers the the second vector/matrix multiplied by A (decision)
 Matrix C = refers to the third vector/matrix subtracted from A*B (threshold)

 Original Code Source:
 http://mytechnotrick.blogspot.com/2013/06/c-program-to-matrix-addition.html
 Modified by Adela

 Merged Code Source from Sally Li:
    Author: Ruiqi (Sally) Li
            Bowdoin College '18
 
        Email: rli@bowdoin.edu
 
        Description:
            Input: 
                define the number of agents, the probability of
                an edge forming, the mean and the standard
                deviation for the influence factor distribution,
                and the number of digits desired for the
                influence factors
            Output: 
                a .csv file containing the influence matrix
                a .dot file for gephi to show network structure 
 ******************************************************************************/

/*****************************************************************************/
/* include files */
#include <iostream> 
#include <vector>
#include <string>
#include <cstdlib> 
#include <map>
#include <random>
#include <ctime>
#include <fstream>
#include <sstream>
#include <algorithm> 
#include <iterator>
#include <stack>
#include <iomanip>

using namespace std;

/***************************************************************************/
/* constants */
// CHANGE ANYTHING HERE IN THIS BOX

#define s 100 //maximum dimensions of the matrix

const int MAX_TIME = 300; // total run time of program, 5 minutes = 300 seconds
const int AGENTS = 10;  // actual dimensions of matrix A
const int MAX_ITERATION = 10; // the most amount of iteration allowed until program stops

const double MEAN = 0.5; // mean of normal distribution
const double SD = 0.5;  // standard deviation

// define the number of agents, edge probability,
// influence factor distributions and the number
// of digit desired for the influence factor
//Sally
const double EDGEP = 0.5;
const int PRECISION = 3;
int POSITIVE = 0; //count of positive 1 in matrix C

// create a csv file
ofstream myfile;

/***************************************************************************/
/* FUNCTIONS */

//Sally
void printM(int n);
double b[s][s];
double c[s][s];

/***************************************************************************/
/* MATRIX CLASS */
class matrix
{

/***************************************************************************/
/* globals variables */
double a[s][s];
int x;
int y;
static int n;
  
/***************************************************************************/
/* functions */  
public:
    void get(); // gets matrix from terminal input
    void put(); // prints matrix

    //Sally
    void generateRand(int matrixType);

    void readIn();  // reads in matrix from csv file
    void initA();   // influence factor matrix
    void initB();   // agent decision at time t matrix
    void initC();   // threshold matrix
    matrix updateB(); // converts output to matrix B
    void updateArray(int n); // adds outcomes to the array    
    void printAC(int n); // prints A and C matrix to determine if equilibrium
    bool operator==(matrix);
    matrix operator+(matrix);
    matrix operator-(matrix);
    matrix operator*(matrix);
    matrix transpose();
};

void matrix:: printAC(int n){
    if(n==0){
        for(int i = 0; i < x; i++){
            for(int j = 0; j < y; j++){
                c[i][j] = a[i][j];
            }  
        }
    }
    else {
        for(int i = 0; i < x; i++){
            c[i][n] = a[i][0];
        }

        myfile.open("InfluenceMatrix.txt");

        myfile << AGENTS << endl;

        // print each influence factor
        for (int i = 0; i < AGENTS+1; i++){
            for (int j = 0; j < AGENTS+1; j++){
                myfile << c[i][j] << " ";
            }
            
            myfile << endl;
        }

        myfile.close();
    }

}

/*****************************************************************************
 Function:  updateArray
 Inputs:    nothing
 Returns:   nothing
 Description: adds outcomes to the array
 *****************************************************************************/
void matrix:: updateArray(int n){
    for(int i = 0; i < x; i++){
        b[i][n] = a[i][0];
    }
}

/*****************************************************************************
 Function:  updateB
 Inputs:    nothing
 Returns:   updated matrix B
 Description: converts answer from calculations to updated matrix B
 *****************************************************************************/
matrix matrix::updateB(){
    matrix r;
    for(int i = 0; i < x; i++) {
        for(int j = 0; j < y; j++) {
            if(a[i][j] < 0){
                r.a[i][j] = -1;
            }
            else if(a[i][j] > 0){
                r.a[i][j] = 1;
            }
            else{
                r.a[i][j] = 0;
            }
        }
    }

    r.x = x;
    r.y = y;

    return r;
}

/*****************************************************************************
 Function:  generateRand
 Inputs:    none
 Returns:   none
 Description: initialize matrix A with a normal distribution and erdos renyi
 *****************************************************************************/
void matrix::generateRand(int matrixType) {
    x = AGENTS;
    
    switch (matrixType) {
        case 1: //A
            y = AGENTS;
            break;
        case 2: //B
            y = 1;
            break;
        case 3: //C
            y = 1;
            break;
        default: 
            break;
    }

    // initialize random variable distribution generator
    default_random_engine generator(time(NULL));

    // input for bernoulli: the probability of an edge forming between two nodes
    // input for normal: mean and standard deviation
    bernoulli_distribution distribution1(EDGEP);
    normal_distribution<double> distribution2(MEAN, SD);

    if (matrixType == 3) {
        for (int i = 0; i < x; i++){
            for (int j = 0; j < y; j++){
                a[i][j] = 0;
            }
        }
        return;
    }

    if (matrixType == 2) {
        for (int i = 0; i < x; i++){
            for (int j = 0; j < y; j++){
                a[i][j] = -1;
                if(a[i][j]==1){
                    POSITIVE++;
                }
            }
        }
        return;
    }

    // go through the entire influence matrix
    for (int i = 0; i < x; i++){
        for (int j = 0; j < y; j++){

            // initialize each entry
            a[i][j] = 0;

            // if it is not a diagonal entry
            // else{
                if (i != j){

                    // true value = there is an edge (with probability edgeP)
                    // false value = there does not exist an edge
                    // execute the codes if there is a edge between the 2 nodes
                    if (distribution1(generator)){

                        // generate an influence factor between -1 and 1
                        // the factors will follow a normal distribution
                        double influence;
                        do {
                            influence = distribution2(generator);
                        } while (influence > 1 || influence < -1);

                        // assign that value to the influence matrix
                        a[i][j] = influence;
                    }
                }
            // }
        }
    }
}

/*****************************************************************************
 Function:  readIn
 Inputs:    
 Returns:   matrix
 Description: reads in content of file and stores as a matrix
 *****************************************************************************/
void matrix::readIn(){

    string fileName;
    ifstream dataFile;
    vector<double> tempMatrix;

    cout << "Enter the number of rows:" << endl;
    cin >> x;

    cout << "Enter the number of columns:" << endl;
    cin >> y;

    // open file
    do {
        cout << "Enter a csv file:" << endl;
        cin >> fileName;
        dataFile.open(fileName.c_str(), ios::in);
    } while (!dataFile.good());

    string line;

    // check that file is not empty
    // and also get rid of the first line
    // that are just data type names
    if (!getline(dataFile, line)) {
        cerr << "See error: The csv file is empty" << endl;
        return;
    }

    // indicate that the file is read in

    //get data of the first line
    stringstream  lineStream(line);
    string cell;
        
    while(getline(lineStream,cell,',')){
        tempMatrix.push_back(stod(cell.c_str()));
    }
    
    //gets data of the rest of the lines
    while(getline(dataFile,line)) {
        stringstream  lineStream(line);
        string cell;
        
        while(getline(lineStream,cell,',')){
            tempMatrix.push_back(stod(cell.c_str())); 
        }
    }

    if(tempMatrix.size() != (x*y)) {
        cout << "These entries do not match the dimensions of this matrix." << endl;
        exit(1);
    }

    int k = 0;
    for(int i = 0; i < x; i++){ 
        for(int j = 0; j < y; j++){ 
            a[i][j] = tempMatrix[k];
            k++;
        }
    }    

    dataFile.close();
}

/*****************************************************************************
 Function:  initA
 Inputs:    nothing
 Returns:   nothing
 Description: initialize matrix A
 *****************************************************************************/
void matrix::initA() {

    x = AGENTS;
    y = AGENTS;

    double tempA[] = {0, 0.9, 0.1, 0.4, 0, 0.1, 0, 0, 0, 0};

    int k = 0;
    for(int i = 0; i < x; i++){ 
        for(int j = 0; j < y; j++){ 
            a[i][j] = tempA[k];
            k++;
        }
    }
}

/*****************************************************************************
 Function:  initB
 Inputs:    nothing
 Returns:   nothing
 Description: initializes matrix B
 *****************************************************************************/
void matrix::initB() {

    x = AGENTS;
    y = 1;

    //CHANGE
    double tempA[] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

    int k = 0;
    for(int i = 0; i < x; i++){ 
        for(int j = 0; j < y; j++){ 
            a[i][j] = tempA[k];
            if(a[i][j] == 1){
                POSITIVE++;
            }
            k++;
        }
    }
}

/*****************************************************************************
 Function:  initC
 Inputs:    nothing
 Returns:   nothing
 Description: initializes matrix C
 *****************************************************************************/
void matrix::initC() {

    x = AGENTS;
    y = 1;

    //CHANGE
    double tempA[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};


    int k = 0;
    for(int i = 0; i < x; i++){ 
        for(int j = 0; j < y; j++){ 
            a[i][j] = tempA[k];
            k++;
        }
    }
}

/*****************************************************************************
 Function:  get
 Inputs:    nothing
 Returns:   nothing
 Description: reads in data of matrix
 *****************************************************************************/
void matrix::get() {
    cout << "Enter the order of Matrix "<<" :\n";
    cin >> x >> y;
    cout << "Enter the Matrix " << " :\n";
    for(int i = 0; i < x; i++) 
        for(int j = 0; j < y; j++) 
            cin >> a[i][j];
}

/*****************************************************************************
 Function:  put
 Inputs:    nothing
 Returns:   nothing 
 Description: prints out matrix
 *****************************************************************************/
void matrix::put() {
    for(int i = 0; i < x; i++) {
        cout << "\n\t";
        myfile << endl; 
        for(int j = 0; j < y; j++) {
            cout << a[i][j] << " ";
            myfile << a[i][j] << ",";
        }
    }

    cout << "\n";
    myfile << endl;
}

/*****************************************************************************
 Function:  ==
 Inputs:    matrix
 Returns:   matrix
 Description: checks to see if two matrices are equal
 *****************************************************************************/
bool matrix::operator==(matrix b) {
    if( (x != b.x) || (y != b.y)) {
        cout << "\n\tMatrix equality is not possible the result is incorrect\n\n";
    }

    for(int i = 0; i < x; i++) {
        for(int j = 0; j < y; j++) {
            if(a[i][j] != b.a[i][j]){
                return false;
            }
        }
    }

    return true;
}

/*****************************************************************************
 Function:  +
 Inputs:    matrix
 Returns:   matrix
 Description: adds two matrices
 *****************************************************************************/
matrix matrix::operator+(matrix b) {
    matrix r;
    if( (x != b.x) || (y != b.y)) {
        cout << "\n\tMatrix Addition is not possible the result is incorrect\n\n";
        r.x = 0;
        r.y = 0;
    }
    else {
        r.x = x;
        r.y = y;
    }

    for(int i = 0; i < x; i++)
        for(int j = 0; j < y; j++)
            r.a[i][j] = a[i][j] + b.a[i][j];

    return r;
}

/*****************************************************************************
 Function:  -
 Inputs:    matrix
 Returns:   matrix
 Description: subtracts two matrices
 *****************************************************************************/
matrix matrix::operator-(matrix b){
    matrix r;
    if( (x != b.x) || (y != b.y) ) {
        cout << "\n\tMatrix subtraction is not possible the result is incorrect\n\n";
        r.x = 0;
        r.y = 0;
    }
    else {
        r.x = x;
        r.y = y;
    }
    for(int i = 0; i < x; i++)
        for(int j = 0; j < y; j++)
            r.a[i][j] = a[i][j] - b.a[i][j];

    return r;
}

/*****************************************************************************
 Function:  *
 Inputs:    matrix
 Returns:   matrix
 Description: multiplies two matrices together
 *****************************************************************************/
matrix matrix::operator*(matrix b){
    matrix r;
    if( (y != b.x) ) {
        cout << "x " << x << " b.y " << b.y << " y "  << y << " b.x "<< b.x << endl;
        cout << "\n\tMatrix Multiplication is not possible the result is incorrect\n\n";
        r.x = 0;
        r.y = 0;
    }
    else {
        r.x = x;
        r.y = b.y;
    }
    
    int i, j, k;
    for(i = 0; i < s; i++)
        for(j = 0; j < s; j++)
            r.a[i][j] = 0;
    for(i = 0; i < x; i++)
        for(j = 0;j < b.y; j++)
            for(k = 0; (k < y) || (k < b.x); k++)
                r.a[i][j] += a[i][k] * b.a[k][j];

    return r;
}

/*****************************************************************************
 Function:  transpose
 Inputs:    nothing
 Returns:   transpose of matrix
 Description: transposes the matrix
 *****************************************************************************/
matrix matrix::transpose()
{
    matrix r;
    for(int i = 0;i < x; i++)
        for(int j = 0;j < y; j++)
            r.a[i][j] = a[j][i];

    r.x = x;
    r.y = y;

    return r;
}

/*****************************************************************************
 Function:  print
 Inputs:    nothing
 Returns:   nothing
 Description: nothing
 *****************************************************************************/
void printM(int n){

    // create a csv file
    // ofstream myfile;

    // turn on if you want to limit the number of
    // digits in the influence factor
//  myfile<<setprecision(precision);

    for(int i=0; i <= n; i++) {
        cout << " "  << i << " ";
        myfile << i << ",";
    }

    cout << endl;
    myfile << endl;
    // print each influence factor
    for (int i = 0; i < AGENTS; i++){
        for (int j = 0; j <= n; j++){
        // int j=0;

            // the comma at the end meanfƒs it is the end of
            // the values for this cell
           // myfile << arr[i][j] << ",";
            if(b[i][j] < 0){
                myfile << b[i][j] << ",";
                cout << b[i][j] << " ";
            }
            else {
                myfile << b[i][j] << ",";
                cout << " " << b[i][j] << " ";
            }
        }

        // goto the next line for the next person
        myfile << endl;
        cout << endl;
    }

    // myfile.close();
}


/*****************************************************************************
 Function:  main
 Inputs:    nothing
 Returns:   matrix answers
 Description:   
            runs the algorithm for doing matrix operations and produces answers
 *****************************************************************************/
int main(){
    matrix a, b, c, d, e;
    int n = 0; // time iteration
    myfile.open(to_string(AGENTS) + " " + to_string(POSITIVE) + " " + to_string(MEAN) + " " + to_string(EDGEP) + ".csv");

    // generate matrix A using distribution
    //CHANGE
    a.generateRand(1);
    cout << "This is your matrix A: Influence factors " << endl;
    myfile << "This is your matrix A: Influence factors " << endl;
    a.put();
    // e = a.transpose();
    // e.printAC(0);
    // generate matrix B using human input
    // if want distribution to be initialized all to -1, use b.generateRand(2)
    //CHANGE
    // b.initB();
    b.generateRand(2);
    cout << "This is your matrix B: Agents decision at time t=" << n << endl;
    myfile << "This is your matrix B: Agents decision at time t=" << n << endl;
    b.put();
    // generate matrix C using human input
    // if want distribution to be initialized all to 0, use c.generateRand(3)
    //CHANGE
    // c.initC();
    c.generateRand(3);
    cout << "This is your matrix C: Threshold " << endl;
    myfile << "This is your matrix C: Threshold " << endl;
    c.put();
    myfile << endl;
    // c.printAC(AGENTS);

    d = (a * b) - c;
    b = d.updateB();
    // b.put();
    b.updateArray(n);
    // printM(n);

    // If entry in matrix c ≥ 0 -> entry in b = -1. If entry in matrix c < 0 -> entry in b = 1.
    matrix previous;
    previous = b;
    clock_t begin = clock();
    clock_t end;
    double elapsed_secs = 0;

    //CHANGE
    // can choose either option of total time, or max iterations
    // while(elapsed_secs < MAX_TIME){
    while(n < MAX_ITERATION){
        n++;
        d = (a * b) - c;

        b = d.updateB();
        cout << "This is your matrix Updated B: Agents decision at time t=" << n << endl;
        // b.put();
        b.updateArray(n);
        // printM(n);
        

        if(previous == b){
            break;
        }
        previous = b;
        end = clock();
        elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    }
    printM(n);
    myfile.close();
        
    return 0;
}

/****** END OF FILE **********************************************************/