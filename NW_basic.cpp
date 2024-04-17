#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
#include <chrono> // Include the <chrono> header for timing
#include <iomanip> // Include for std::setw

using namespace std;
using namespace std::chrono; // Use the std::chrono namespace for timing

// Function to compute the Needleman-Wunsch score matrix and trace back for alignment
tuple<vector<vector<int> >, string, string> NW_algorithm(const string& seq1, const string& seq2, int match, int mismatch, int gap) {
    auto start = high_resolution_clock::now(); // Start timing

    int len1 = seq1.length();
    int len2 = seq2.length();
    vector<vector<int> > score(len1 + 1, vector<int>(len2 + 1));
    vector<vector<pair<int, int> > > trace(len1 + 1, vector<pair<int, int> >(len2 + 1));

    // Initialization
    for (int i = 0; i <= len1; i++) {
        score[i][0] = i * gap;
        if (i > 0) trace[i][0] = make_pair(i - 1, 0);
    }
    for (int j = 0; j <= len2; j++) {
        score[0][j] = j * gap;
        if (j > 0) trace[0][j] = make_pair(0, j - 1);
    }

    // Filling the score matrix and trace
    for (int i = 1; i <= len1; i++) {
        for (int j = 1; j <= len2; j++) {
            int matchScore = score[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match : mismatch);
            int deleteScore = score[i - 1][j] + gap;
            int insertScore = score[i][j - 1] + gap;
            score[i][j] = max(max(matchScore, deleteScore), insertScore);
            if (score[i][j] == matchScore) {
                trace[i][j] = make_pair(i - 1, j - 1);
            } else if (score[i][j] == deleteScore) {
                trace[i][j] = make_pair(i - 1, j);
            } else {
                trace[i][j] = make_pair(i, j - 1);
            }
        }
    }

    // Trace back for alignment
    string align1 = "", align2 = "";
    for (int i = len1, j = len2; i > 0 || j > 0;) {
        if (trace[i][j] == make_pair(i - 1, j - 1)) {
            align1 = seq1[i - 1] + align1;
            align2 = seq2[j - 1] + align2;
            i--;
            j--;
        } else if (trace[i][j] == make_pair(i - 1, j)) {
            align1 = seq1[i - 1] + align1;
            align2 = "-" + align2;
            i--;
        } else {
            align1 = "-" + align1;
            align2 = seq2[j - 1] + align2;
            j--;
        }
    }

    auto stop = high_resolution_clock::now(); // Stop timing
    auto duration = duration_cast<microseconds>(stop - start); // Calculate duration
    cout << "Time taken by function: " << duration.count() << " microseconds" << endl;

    return make_tuple(score, align1, align2); // Use make_tuple to create a tuple object
}

// Function to print the score matrix and alignment with better alignment across rows and columns
void printMatrixAndAlignment(const vector<vector<int> >& matrix, const string& align1, const string& align2) {
    int max_width = 0;
    for (const auto& row : matrix) {
        for (const auto& elem : row) {
            int elem_width = to_string(elem).length() + 1; // +1 for space
            max_width = max(max_width, elem_width);
        }
    }

    for (const auto& row : matrix) {
        for (const auto& elem : row) {
            cout << setw(max_width) << elem;
        }
        cout << endl;
    }
    cout << endl;
    cout << "Alignment:" << endl;
    cout << align1 << endl;
    cout << align2 << endl;
}

int main() {
    string seq1 = "TAATTGTATCCCAT";
    string seq2 = "GACGACGACGATCA";
    int match = 2;
    int mismatch = -1;
    int gap = -2;

    // print the seq1 and seq2
    cout << "Sequence 1: " << seq1 << endl;
    cout << "Sequence 2: " << seq2 << endl;
    cout << endl;

    tuple<vector<vector<int> >, string, string> result = NW_algorithm(seq1, seq2, match, mismatch, gap); // Correctly declare the tuple type
    vector<vector<int> > scoreMatrix = get<0>(result);
    string align1 = get<1>(result);
    string align2 = get<2>(result);
    cout << "Needleman-Wunsch Score Matrix and Alignment:" << endl;
    printMatrixAndAlignment(scoreMatrix, align1, align2);

    return 0;
}
