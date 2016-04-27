#include <iostream>
#include <map>

template<typename TMap>
void func(TMap) {
    TMap m;
}

int main(int argc, char *argv[])
{
    typedef std::map<int, int> TMap;
    TMap m;
    func(TMap());
    return 0;
}

