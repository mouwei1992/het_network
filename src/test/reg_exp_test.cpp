#include<iostream>
#include<regex>
#include<string>

using namespace std;

int main(){
  string test_string = "11 10  00 0 1";
  regex num("([0-9]+)");

  sregex_iterator currentMatch(test_string.begin(),test_string.end(), num);
  sregex_iterator lastMatch;
  smatch matches;
  int n_matches = 0;

  while(currentMatch != lastMatch){
       std::smatch match = *currentMatch;
       std::cout << match.str() << "\n";
       currentMatch++;
       n_matches++;
   }

   cout << n_matches << endl;

  return 0;
}
