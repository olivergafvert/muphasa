#ifndef MPH_GRADE_INCLUDED
#define MPH_GRADE_INCLUDED

#include "utils.h"
#include <vector>
#include <queue>
#include <iostream>
#include <set>

// TODO: define object monomialOrder, it implements leq

class grade_t : public std::vector<index_t> {
 public:
    // TODO: alternative use Array which implements some lexicographical comparison
    // https://en.cppreference.com/w/cpp/container/array/operator_cmp
    
  using std::vector<index_t>::vector;
  
  grade_t(std::vector<index_t>& v) : std::vector<index_t>(v){};
  
  grade_t(const grade_t& v) : std::vector<index_t>(v){};
    
  bool leq_poset(grade_t& other_grade); /* Poset order */
    
  bool leq_poset_m(grade_t& other_grade);
    
  bool lt_poset(grade_t& other_grade); /* Poset order */
  
  bool lt_colex(const grade_t& other_grade) const; /* Colexicographical order */
  
  bool operator==(grade_t& other_grade) const;
    
  bool operator<(grade_t& other_grade) const; /* Lexicographical order */
  
  grade_t join(grade_t& other_grade);
    
  grade_t m_ji(grade_t& other_grade);
    
  void print();
    
    std::string to_string();
};

inline bool grade_t::leq_poset(grade_t& other_grade)
{
/* Poset order */
  for(size_t i=0; i<other_grade.size(); i++){
    if(this->at(i) > other_grade[i]){
      return false;
    }
  }
  return true;
}

inline bool grade_t::leq_poset_m(grade_t& other_grade)
{
/* Poset order */
  for(size_t i=0; i<other_grade.size(); i++){
    if(this->at(i) > other_grade[i]){
      return false;
    }
  }
  return true;
}

inline bool grade_t::lt_poset(grade_t& other_grade)
{
    /* Poset order */
    bool is_lt = false;
    for(size_t i=0; i<other_grade.size(); i++){
        if(this->at(i) > other_grade[i]){
            return false;
        }else if(this->at(i) < other_grade[i]){
            is_lt = true;
        }
    }
    return is_lt;
}

inline bool grade_t::lt_colex(const grade_t& other_grade) const
{
/* Colexicographical order */
    if(this->size()!=other_grade.size()){
        std::cout<<"error";
    }
    assert(this->size()==other_grade.size());
    size_t i=other_grade.size()-1;
    while (i>0 && other_grade[i]==this->at(i)){
        i--;
    }
    return this->at(i)<other_grade[i];
}



inline bool grade_t::operator==(grade_t& other_grade) const
{
    size_t s = other_grade.size();
    for(size_t i=0; i<s; i++){
        if(other_grade[i] != this->at(i)){
            return false;
        }
    }
    return true;
}

inline bool grade_t::operator<(grade_t& other_grade) const
{
    /* Lexicographical order */
    size_t i=0;
    while(i<other_grade.size() && this->at(i) == other_grade[i]){
        i++;
    }
    if(i==other_grade.size()){
        return false;
    }
    return this->at(i) < other_grade[i];
}

inline grade_t grade_t::join(grade_t& other_grade)
{
    grade_t grade(other_grade);
    for(size_t i=0; i<other_grade.size(); i++){
        if(this->at(i) > grade[i]){
            grade[i] = this->at(i);
        }
    }
    return grade;
}

inline grade_t grade_t::m_ji(grade_t& other_grade)
{
    grade_t grade(other_grade);
    for(size_t i=0; i<other_grade.size(); i++){
        if(this->at(i) < grade[i]){
            grade[i] = 0;
        }else{
            grade[i] = this->at(i) - grade[i];
        }
    }
    return grade;
}

inline void grade_t::print(){
  if(this->size()>0){
    std::cout << "(" << this->at(0);
    for(size_t i=1; i<this->size(); i++){
      std::cout << ", " << this->at(i);
    }
    std::cout << ")\n";
  } else{
    std::cout << "()\n";
  }
  std::cout << std::flush;
}

inline std::string grade_t::to_string(){
    std::string output;
    if(this->size()>0){
        output += "(";
        output += std::to_string(this->at(0));
        for(size_t i=1; i<this->size(); i++){
            output += ", ";
            output += std::to_string(this->at(i));
        }
        output += ")";
    } else{
        output += "()";
    }
    return output;
}

class Latt_iterator{
private:
    std::vector<grade_t> grades;
    std::vector<grade_t> indices;
    
    int curr_i=0;
    int indices_index = 0;
    
    void gen_rec(int i, int r, std::vector<grade_t> _grades){
        std::vector<grade_t> new_grades;
        new_grades.reserve(i);
        for(int j=0; j<=i; j++){
            new_grades.push_back(_grades[i].join(_grades[j]));
        }
        sort(new_grades.begin(), new_grades.end());
        if(r==1){
            if(indices.size() == 0 || new_grades[0] != indices.back()){
                indices.push_back(new_grades[0]);
            }
            for(int j=1; j<new_grades.size(); j++){
                if(new_grades[j] != indices.back()){
                    indices.push_back(new_grades[j]);
                }
            }
        }else{
            for(int j=0; j<new_grades.size(); j++){
                gen_rec(j, r-1, new_grades);
            }
        }
    }
public:
    
    Latt_iterator(std::vector<grade_t>& _grades){
        grades = _grades;
        sort(grades.begin(), grades.end()); // Sort grades in lex order
    }
    
    grade_t& next(){
        if(indices_index == indices.size()){
            indices.clear();
            indices_index = 0;
            gen_rec(curr_i, (int)grades[0].size()-1, grades);
            curr_i++;
        }
        grade_t& ret = indices[indices_index];
        indices_index++;
        return ret;
    }
    
    bool has_next(){
        return !(curr_i == grades.size() && indices_index == indices.size());
    }
};

class Iterator_lex{
public:
    /** Lexicographical iterator over a base set. */
    std::vector<std::vector<index_t>> base_set; // Base set to iterate over
    grade_t indices; // Index vector that specifies current element
    grade_t values; // Values of the current index
    Iterator_lex(const std::vector<std::vector<index_t>>& base_set){
        this->base_set = base_set;
        
        // Initialize the index vector
        for(size_t i=0; i<base_set.size(); i++){
            indices.push_back(0);
            values.push_back(base_set[i][0]);
        }
        if(base_set.size() > 0){
            indices[base_set.size()-1] = -1; //Initial state needs to be negative
        }
    }
    
    bool has_next();
    
    grade_t& next();
    
    void step();
};

inline bool Iterator_lex::has_next(){
  for(size_t i=0; i<indices.size(); i++){
    if(indices[i] == -1 || indices[i]<base_set[i].size()-1){
      return true;
    }
  }
  return false;
}

inline void Iterator_lex::step(){
    index_t k = indices.size()-1;
    while(k>0 && indices[k] == base_set[k].size()-1){
        indices[k] = 0;
        values[k] = base_set[k][0];
        k--;
    }
    indices[k]++;
    values[k] = base_set[k][indices[k]];
}
    
inline grade_t& Iterator_lex::next(){
  index_t k = indices.size()-1;
  while(k>0 && indices[k] == base_set[k].size()-1){
    indices[k] = 0;
    values[k] = base_set[k][0];
    k--;
  }
  indices[k]++;
  values[k] = base_set[k][indices[k]];
  return values;
}

template <typename Column> std::vector<std::vector<index_t>> get_grade_base_set(const std::vector<Column>& columns){
    /** Returns a list of lists where the i-th list contains the unique values of the grades along the i-th coordinate axis.
     
     Arguments:
     columns {std::vector<Column>} -- a list of graded columns.
     
     Returns:
     std::vector<std::vector<index_t>> -- a list of lists where the i-th list contains the unique values of the grades along the i-th coordinate axis.
     */
    std::vector<std::vector<index_t>> index_value_lists;
    if(columns.size() == 0){
        return index_value_lists;
    }
    
    std::vector<std::set<index_t>> index_values;
    for(size_t i=0; i<columns[0].grade.size(); i++){
        index_values.push_back(std::set<index_t>());
    }
    
    for (size_t index_column = 0; index_column < columns.size(); index_column++){
        /* Add indices to be visited to the vector 'indices' */
        for(size_t index_grade=0; index_grade< columns[index_column].grade.size(); index_grade++){
            index_values[index_grade].insert(columns[index_column].grade[index_grade]);
        }
    }
    
    for(size_t i=0; i<index_values.size(); i++){
        index_value_lists.push_back(std::vector<index_t>(index_values[i].begin(), index_values[i].end()));
        sort(index_value_lists[i].begin(), index_value_lists[i].end());
    }
    
    return index_value_lists;
}

std::vector<std::vector<grade_t>> get_grade_slices(const std::vector<std::vector<index_t>>& base_set){
    std::vector<std::vector<grade_t>> grade_slice_list;
    Iterator_lex iterator(base_set);
    while(iterator.has_next()){
        grade_t v = iterator.next();
        index_t sum = 0;
        for(auto& val : v){
            sum += val;
        }
        while(grade_slice_list.size() <= sum){
            grade_slice_list.push_back(std::vector<grade_t>());
        }
        grade_slice_list[sum].push_back(v);
    }
    return grade_slice_list;
}


void helper_grade_slice(grade_t& indices, size_t sum, size_t target_sum, size_t index, std::vector<grade_t*>& grades, std::vector<std::vector<index_t>>& base_set, size_t n_parameters){
    if(index+1==n_parameters){
        if(base_set[index].size()+sum > target_sum){
            std::vector<index_t> grade;
            for(size_t i=0; i<n_parameters-1; i++){
                grade.push_back(base_set[i][indices[i]]);
            }
            grade.push_back(base_set[index][target_sum-sum]);
            grades.push_back(new grade_t(grade));
        }
    }else{
        for(size_t i=0; i<base_set[index].size(); i++){
            if(sum+i <= target_sum){
                indices[index] = i;
                helper_grade_slice(indices, sum+i, target_sum, index+1, grades, base_set, n_parameters);
            }else{
                break;
            }
        }
    }
}

std::vector<grade_t*> get_grade_slice_rec(std::vector<std::vector<index_t>>& base_set, int val){
    std::vector<grade_t*> grade_slice_list;
    grade_t indices;
    for(size_t i=0; i<base_set.size(); i++){
        indices.push_back(0);
    }
    helper_grade_slice(indices, 0, val, 0, grade_slice_list, base_set, base_set.size());
    return grade_slice_list;
}

std::vector<grade_t> get_grade_slice(const std::vector<std::vector<index_t>>& base_set, int val){
    std::vector<grade_t> grade_slice_list;
    Iterator_lex iterator(base_set);
    while(iterator.has_next()){
        iterator.step();
        index_t sum = 0;
        for(auto& val : iterator.indices){
            sum += val;
        }
        if(sum==val){
            grade_slice_list.push_back(grade_t(iterator.values));
        }
    }
    return grade_slice_list;
}


struct GradeHasher
{
    std::size_t operator()(const grade_t& k) const
    {
        return boost::hash_range(k.begin(), k.end());
    }
};

#endif // MPH_GRADE_INCLUDED


