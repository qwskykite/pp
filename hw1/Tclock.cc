#include<ctime>
#include<string>
#include<cstdio>
using namespace std;
class Tclock{
    public:
        Tclock(string title);
        void AddDura();
        void Reset();
        void Print();
        
    private:
        clock_t start;
        double totalTime;
        string title;
};
Tclock::Tclock(string title){
    this->title = title;
    Reset();
    start = clock();
    totalTime = 0;
}
void Tclock::AddDura(){
    totalTime += double(clock() - start) / CLOCKS_PER_SEC;
    start = clock();
}
void Tclock::Reset(){
    start = clock();
    totalTime = 0;
}
void Tclock::Print(){
    printf("%s: %lfs\n", title.c_str(), totalTime);
}
int main(){
    int k = 0;
    Tclock ft = Tclock("fxx");
    for(int i = 0; i < 1e9; i++)
        k = k + 1;
    ft.AddDura();
    ft.Print();
        for(int i = 0; i < 1e9; i++)
        k = k + 1;
    ft.AddDura();
    ft.Print();
    return 0;
}