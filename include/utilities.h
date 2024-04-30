#ifndef UTILITIES_INCLUDED
#define UTILITIES_INCLUDED

#include <cassert>
#include <vector>

using std::vector;

//Get an element index in its vector from its pointer
template <typename T>
int index(const T* const& at, const vector<T>& vec){
    assert((at >= &vec.front() && at <= &vec.back()));
    return static_cast<int>(at - &vec.front());
};

//Get an element index in its vector from its iterator
template <typename T>
int index(const typename vector<T>::iterator& at, const vector<T>& vec){
    assert(at >= vec.begin() && at < vec.end());
    return static_cast<int>(at - vec.begin());
};

//Add an element at the end of vector, return a pointer of the appended
template <typename T>
T* append(const T& t_new, vector<T>& vec, T* const& anchor){
    vec.emplace_back(t_new);
    assert(anchor == &vec.front());
    return &vec.back(); 
};

//Erase an element in the vector, move the last element to the erased place, return a pointer of the erased
template <typename T> 
T* erase(T* const& t_ptr, vector<T>& vec){
    assert(index(t_ptr, vec) != -1 && t_ptr != nullptr);
    *t_ptr = vec.back();
    vec.erase(vec.end()-1);
    return t_ptr;
};

//Return sign of variable
// template <typename T> inline constexpr
// int sgn(T val){
//     return (T(0) < val) - (val < T(0));
// };



#endif