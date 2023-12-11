#include <iostream>
#include <vector>
#include <climits>
#include <algorithm>
#include <numeric>

using namespace std;

void inputMatrix(vector<vector<int>>& matrix, int rows, int cols) {
    cout << "Enter the matrix:\n";
    for (int i = 0; i < rows; ++i) {
        matrix.push_back(vector<int>());
        for (int j = 0; j < cols; ++j) {
            int element;
            cout << "Enter an element [" << i + 1 << "][" << j + 1 << "]: ";
                        cin >> element;
            matrix[i].push_back(element);
        }
    }
}

void printMatrix(const vector<vector<int>>& matrix) {
    for (const auto& row : matrix) {
        for (int element : row) {
            cout << element << "\t";
        }
        cout << endl;
    }
}

void decisionCriteria(const vector<vector<int>>& matrix) {

    cout << endl;

    // MaxiMax
    vector<int> maxiMax;
    for(const auto& row : matrix){
        maxiMax.push_back(*max_element(row.begin(), row.end()));
    }
    int resMaxiMax = 0;
    int indexMaxiMax = 0;

    for(int i = 0; i<maxiMax.size(); i++){
        if(maxiMax[i] > resMaxiMax){
            resMaxiMax = maxiMax[i];
            indexMaxiMax = i;
        }
    }


    // Vald
    vector<int> vald;
    for (const auto& row : matrix) {
        vald.push_back(*min_element(row.begin(), row.end()));
    }
    // Find max
    int res = INT_MIN;
    int index = 0;
    for(int i = 0; i<vald.size(); i++){
        if(vald[i] > res){
            res = vald[i];
            index = i;
        }
    }

    // Sevich
    vector<int> sevich;
    for (int j = 0; j < matrix[0].size(); ++j) {

        // Max from columns
        int maxVal = matrix[0][j];
        for (int i = 1; i < matrix.size(); ++i) {
            maxVal = max(maxVal, matrix[i][j]);
        }
        sevich.push_back(maxVal);
    }
    // Result for Sevich

    // init new matrix
    vector<vector<int>> table(matrix.size(),vector<int>(matrix[0].size(),0));

    // Put matrix elements to new one
    for(int i = 0; i<matrix.size(); i++){
        for(int j = 0; j<matrix[i].size(); j++){
            table[i][j] = matrix[i][j];
        }
    }

    // Calc. differences
    for(int i = 0; i<matrix.size(); i++){
        for(int j = 0; j<matrix[i].size(); j++){
            table[i][j] -= sevich[j];
            table[i][j] = abs(table[i][j]);
        }
    }

    // Find MinMax
    vector<int> resSevich;
    for (const auto& row : table) {
        resSevich.push_back(*max_element(row.begin(), row.end()));
    }
    int resSv = INT_MAX;
    int indexSv = 0;
    for(int i = 0; i<resSevich.size(); i++){
        if(resSevich[i] < resSv){
            resSv = resSevich[i];
            indexSv = i;
        }
    }

    // Gurvitz
    double alpha;
    cout << "Enter the 'alpha' for Gurvitz : ";
    cin >> alpha;
    vector<double> gurvitz;
    for (int i = 0; i < matrix.size(); ++i) {
        double maxVal = *max_element(matrix[i].begin(), matrix[i].end());
        double minVal = *min_element(matrix[i].begin(), matrix[i].end());
        gurvitz.push_back((1 - alpha) * maxVal + alpha * minVal);
    }

    // Find the max. el.
    int resGur = INT_MIN;
    int indexGur = 0;
    for(int i = 0; i<gurvitz.size(); i++){
        if(gurvitz[i] > resGur){
            resGur = gurvitz[i];
            indexGur = i;
        }
    }


    // Laplace
    vector<double> laplace;
    for (const auto& row : matrix) {

        // Sum rows and divide by row.size
        double average = accumulate(row.begin(), row.end(), 0.0) / row.size();
        laplace.push_back(average);
    }

    // Find max
    int resLap = INT_MIN;
    int indexLap = 0;
    for(int i = 0; i<laplace.size(); i++){
        if(laplace[i] > resLap){
            resLap = laplace[i];
            indexLap = i;
        }
    }

    // Bayes
    vector<double> bayes;
    int resBayes = INT_MIN;
    int indexBayes = 0;
    double check = 0;

    cout << "Enter the  probabilities for Bayes:\n";
    vector<double> probabilities(matrix[0].size());

    for (int j = 0; j < matrix[0].size(); ++j) {
        cout << "Probability of " << j + 1 << ": ";
        cin >> probabilities[j];
        check += probabilities[j];
    }
    if(check != 1){
        cout << "\nProbabilities were entered incorrectly" << endl;
    }
    else{

    for (int i = 0; i < matrix.size(); ++i) {
        double sum = 0;
        for (int j = 0; j < matrix[0].size(); ++j) {

            // (el. of the row)*probability
            sum += matrix[i][j] * probabilities[j];
        }
        bayes.push_back(sum);
    }

    // Find max
    for(int i = 0; i<bayes.size(); i++){
        if(bayes[i] > resBayes){
            resBayes = bayes[i];
            indexBayes = i;
        }
    }
    }

    // Hodges-Lehmann

    double coef = 0;
    cout << "Enter the coefficient for Hodges-Lehmann ";
    cin >> coef;

    vector<double> hodgLem;

    // Mat. elem * probabilities
    vector<vector<double>> hodgLemMat(matrix.size(), vector<double>(matrix[0].size(),0));

    for(int i = 0; i<matrix.size(); i++){
        for(int j=0; j<matrix[i].size(); j++){
            hodgLemMat[i][j] = matrix[i][j]*probabilities[j];
        }
    }

    // Sum of each row in ( elem * prob ) matrix
    vector<double> sumRow;

    for(int i = 0; i<hodgLemMat.size(); i++){
        double sum = 0;
        for(int j=0; j<hodgLemMat[i].size(); j++){
            sum += hodgLemMat[i][j];
        }
        sumRow.push_back(sum/hodgLemMat[i].size());
    }

    int k = 0;
    for(const auto& row : matrix){

        // formula
        double res = (1-coef)*(*min_element(row.begin(), row.end())) + coef*sumRow[k];
        k++;
        hodgLem.push_back(res);
    }

    // Max from hodgLem

    int resHodgLem = INT_MIN;
    int indexHodgLem = 0;
    for(int i = 0; i<hodgLem.size(); i++){
        if(hodgLem[i] > resHodgLem){
            resHodgLem = hodgLem[i];
            indexHodgLem = i;
        }
    }


    // Bayes - Laplace
    vector<double> BayesLaplace;
    for(int i = 0; i<matrix.size(); i++){
        double sum = 0;
        for(int j = 0; j<matrix[i].size(); j++){

            // (el. of the row)*probability
            sum += matrix[i][j]*probabilities[j];
        }
        // Divide by row.size
        BayesLaplace.push_back(sum/matrix[i].size());
    }

    // Find max
    int resBayLap = INT_MIN;
    int indexBayLap = 0;
    for(int i = 0; i<BayesLaplace.size(); i++){
        if(BayesLaplace[i] > resBayLap){
            resBayLap = BayesLaplace[i];
            indexBayLap = i;
        }
    }



    // Results

    cout << "MaxiMax: ";
    cout << "Row number is "<< indexMaxiMax+1 << endl;

    cout << endl;

    cout << "Vald: ";
    cout << "Row number is "<< index+1 << endl;

    cout << endl;

    cout << "Sevidj: ";
    cout << "Row number is "<< indexSv+1 << endl;

    cout << endl;

    cout << "Gurvitz (alpha = " << alpha << "): ";
    cout << "Gurvitz: ";
    cout << "Row number is "<< indexGur+1 << endl;

    cout << endl;

    cout << "Laplas: ";
    cout << "Row number is "<< indexLap+1 << endl;

    cout << endl;

    if(check == 1){
        cout << "Bayes: ";
        cout << "Row number is "<< indexBayes+1 << endl;
        cout << endl;

        cout << "Hodges-Lehmann: ";
        cout << "Row number is "<< indexHodgLem+1 << endl;
        cout << endl;

        cout << "Bayes - Laplace: ";
        cout << "Row number is "<< indexBayLap+1 << endl;
        cout << endl;
    }


}

int main() {
    int rows, cols;
    cout << "Enter the number of strategy for player 1: ";
    cin >> rows;
    cout << "Enter the number of strategy for player 2: ";
    cin >> cols;

    // init for test
    vector<vector<int>> payoffMatrix;

    inputMatrix(payoffMatrix, rows, cols);

    cout << "\nThe matrix:\n";
    printMatrix(payoffMatrix);

    decisionCriteria(payoffMatrix);

    return 0;
}
