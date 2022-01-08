#ifndef MPH_COLUMN_INCLUDED
#define MPH_COLUMN_INCLUDED

#include <queue>
#include <iostream>
#include <utility>

#include "utils.h"

// TODO: This should be of type template (implemenet inside class Column)

class column_entry_t : public std::pair<value_t, index_t> {
    /** Column entry type for sparse column representation. Encodes a tuple (value, position). **/
 public:
  column_entry_t(value_t a, index_t b) : std::pair<value_t, index_t>(a,b) {}
  
  value_t get_value() const { return this->first; }
  index_t get_index() const { return this->second; }
};



template <typename Entry>
class smaller_index {
 public:
    /** Comparator for ordering indices in the columns. */
    bool operator()(const Entry& a, const Entry& b) {
        return a.get_index() < b.get_index();
    }
};

// TODO: Eventually, think about https://en.wikipedia.org/wiki/Skip_list for the column.

class DenseColumn : public std::vector<bool>{
public:
    int local = 0;
    bool modified = true;
    column_entry_t cached_pivot = column_entry_t(0, -1);
    
    DenseColumn() : std::vector<bool>(){};
    
    DenseColumn(int n) : std::vector<bool>(n, false){};
    
    DenseColumn(const DenseColumn& column) : std::vector<bool>(column) {
        cached_pivot = column.cached_pivot;
    }
    
    column_entry_t pop_pivot();
    
    column_entry_t get_pivot();
    
    void push(const column_entry_t& entry);
    
    void plus(DenseColumn& other_column);
    
    void plus_popped(DenseColumn& other_column);
    
    void print();
    
    column_entry_t get_cached_pivot();
    
    void cache_pivot();
    
    void refresh(){return;};
    
    bool empty(){
        if(this->modified){
            return get_pivot().get_index() == -1;
        }
        return cached_pivot.get_index() == -1;
    }
};

inline column_entry_t DenseColumn::pop_pivot(){
    if(this->modified){
        for(int i=(int)this->size()-1; i>=0 ;i--){
            if(this->at(i)){
                this->cached_pivot = column_entry_t(1, i);
                this->modified = false;
                return this->cached_pivot;
            }
        }
        this->cached_pivot = column_entry_t(0, -1);
        this->modified = false;
        return this->cached_pivot;
    }else{
        this->modified = true;
        if(cached_pivot.get_index() != -1){
            this->at(cached_pivot.get_index()) = false;
        }
        return cached_pivot;
    }
}

inline column_entry_t DenseColumn::get_pivot()
{
    if(this->modified){
        column_entry_t result = pop_pivot();
        if (result.get_index() != -1) this->push(result);
        this->cached_pivot = result;
        this->modified = false;
        return result;
    }else{
        return this->cached_pivot;
    }
}

inline void DenseColumn::push(const column_entry_t& entry){
    this->modified = true;
    this->at(entry.get_index()) = entry.get_value()==0 ? false : true;
}

inline void DenseColumn::plus(DenseColumn& other_column){
    this->modified = true;
    for(size_t i=0; i<other_column.size(); i++){
        this->at(i) = this->at(i)!=other_column[i];
    }
}

inline void DenseColumn::plus_popped(DenseColumn& other_column){
    return this->plus(other_column);
}

inline column_entry_t DenseColumn::get_cached_pivot(){ return cached_pivot;}

inline void DenseColumn::cache_pivot(){
    cached_pivot = this->get_pivot();
}

class Column : public std::vector<index_t>{
private:
    size_t inserts_since_refresh = 0; // Threshold of the size of the column for when it will be refreshed.
    
    index_t pop_pivot_index();
public:
    int local = 0;
    bool modified = true;
    column_entry_t cached_pivot = column_entry_t(0, -1);
    
    Column() : std::vector<index_t>() { }
    Column(int n) : std::vector<index_t>() { }
    
    Column(const Column& column) : std::vector<index_t>(column) {
        cached_pivot = column.cached_pivot;
    }
    
    column_entry_t pop_pivot();
  
    index_t get_pivot_index();
    
    column_entry_t get_pivot();
    
    void push(const column_entry_t& entry);
    
    // TODO: Is there a way of merge heaps in CPP in O(n)? [heapification]
    void plus(Column& other_column);
    
    void plus_popped(Column& other_column);
  
    void print();
    
    column_entry_t get_cached_pivot();
    
    void cache_pivot();
    
    void refresh();
};


inline column_entry_t Column::pop_pivot()
{
  this->modified = true;
  if (this->empty()) return column_entry_t(0, -1);
  
  auto pivot = this->front();
  std::pop_heap(this->begin(), this->end());
  this->pop_back();
  while (!this->empty() && this->front() == pivot) {
      std::pop_heap(this->begin(), this->end());
      this->pop_back();
    if (this->empty())
      return column_entry_t(0, -1);
    else {
      pivot = this->front();
        std::pop_heap(this->begin(), this->end());
        this->pop_back();
    }
  }
  return column_entry_t(1, pivot);
}

inline index_t Column::pop_pivot_index()
{
    this->modified = true;
    if (this->empty()) return -1;
    
    auto pivot = this->front();
    std::pop_heap(this->begin(), this->end());
    this->pop_back();
    while (!this->empty() && this->front() == pivot) {
        std::pop_heap(this->begin(), this->end());
        this->pop_back();
        if (this->empty())
            return -1;
        else {
            pivot = this->front();
            std::pop_heap(this->begin(), this->end());
            this->pop_back();
        }
    }
    return pivot;
}

inline column_entry_t Column::get_pivot()
{
  if(this->modified){
      column_entry_t result = pop_pivot();
      if (result.get_index() != -1) this->push(result);
      this->modified = false;
      cached_pivot = result;
      return result;
  }else{
      if (this->empty())
          return column_entry_t(0, -1);
      return cached_pivot;
  }
}

inline index_t Column::get_pivot_index()
{
    return this->get_pivot().get_index();
}

inline void Column::push(const column_entry_t& entry){
    this->modified=true;
    this->push_back(entry.get_index());
    std::push_heap( this->begin( ), this->end( ) );
    return;
}


inline void Column::plus(Column& other_column){
    if(other_column.size()==0){
        return;
    }
    for(index_t entry : other_column){
        this->push_back(entry);
        std::push_heap( this->begin( ), this->end( ) );
    }
    this->inserts_since_refresh += other_column.size();
    this->modified = true;
    if(this->size() < 2*this->inserts_since_refresh){
        this->refresh();
    }
    return;
}

inline void Column::plus_popped(Column& other_column){
    if(other_column.empty()){
        return;
    }
    for(int i=1; i<other_column.size(); i++){
        this->push_back(other_column[i]);
        std::push_heap( this->begin( ), this->end( ) );
    }
    this->inserts_since_refresh += other_column.size()-1;
    this->modified = true;
    if(this->size() < 2*this->inserts_since_refresh){
        this->refresh();
    }
    return;
}


inline void Column::print(){
  for(index_t entry : *this){
      std::cout << entry << ";";
  }
  std::cout << std::endl;
  return;
}

inline column_entry_t Column::get_cached_pivot(){ return cached_pivot;}

inline void Column::cache_pivot(){
    cached_pivot = this->get_pivot();
}

inline void Column::refresh(){
    std::vector<index_t> tmp;
    tmp.reserve(this->size());
    index_t max_element = this->pop_pivot_index();
    while(max_element != -1){
        tmp.push_back(max_element);
        max_element = this->pop_pivot_index();
    }
    std::reverse(tmp.begin(), tmp.end());
    std::make_heap(tmp.begin(), tmp.end());
    this->swap(tmp);
    this->inserts_since_refresh=0;
    return;
}





class VectorColumn : public std::vector<index_t>{
private:
    size_t inserts_since_refresh = 0; // Threshold of the size of the column for when it will be refreshed.
public:
    int local = 0;
    column_entry_t cached_pivot = column_entry_t(0, -1);
    
    VectorColumn() : std::vector<index_t>() { }
    VectorColumn(int n) : std::vector<index_t>() { }
    
    VectorColumn(const VectorColumn& column) : std::vector<index_t>(column) {
        cached_pivot = column.cached_pivot;
    }
    
    VectorColumn(const std::vector<index_t>& column) : std::vector<index_t>(column) {
    }
    
    column_entry_t pop_pivot();
  
    index_t get_pivot_index();
    
    column_entry_t get_pivot();
    
    void push(const column_entry_t& entry);
    
    void plus(VectorColumn& other_column);
    
    void plus_popped(VectorColumn& other_column);
  
    void print();
    
    column_entry_t get_cached_pivot();
    
    void cache_pivot();
    
    void refresh();
};


inline column_entry_t VectorColumn::pop_pivot()
{
    column_entry_t pivot = this->get_pivot();
    if(pivot.get_index() != -1){
        this->pop_back();
    }
    return pivot;
}

inline index_t VectorColumn::get_pivot_index()
{
    return this->empty() ? -1 : this->back();
}

inline column_entry_t VectorColumn::get_pivot()
{
    return this->empty() ? column_entry_t(1, -1) : column_entry_t(1, this->back());
}

inline void VectorColumn::push(const column_entry_t& entry){
    VectorColumn v;
    v.push_back(entry.get_index());
    this->plus(v);
    return;
}

inline void VectorColumn::plus(VectorColumn& other_column){
    VectorColumn temp_col;
    temp_col.local = this->local;
    size_t new_size = this->size() + other_column.size();
    temp_col.resize(new_size);
    std::vector<index_t>::iterator col_end = std::set_symmetric_difference( this->begin(), this->end(), other_column.begin(), other_column.end(), temp_col.begin() );
    temp_col.erase(col_end, temp_col.end());
    this->swap(temp_col);
    return;
}

inline void VectorColumn::plus_popped(VectorColumn& other_column){
    if(other_column.empty()){
        return;
    }
    VectorColumn temp_col;
    temp_col.local = this->local;
    size_t new_size = this->size() + other_column.size()-1;
    temp_col.resize(new_size);
    std::vector<index_t>::iterator col_end = std::set_symmetric_difference( this->begin(), this->end(), other_column.begin(), other_column.end()-1, temp_col.begin() );
    temp_col.erase(col_end, temp_col.end());
    this->swap(temp_col);
    return;
}


inline void VectorColumn::print(){
  for(index_t entry : *this){
      std::cout << entry << ";";
  }
  std::cout << std::endl;
  return;
}

inline column_entry_t VectorColumn::get_cached_pivot(){ return this->get_pivot();}

inline void VectorColumn::cache_pivot(){
    cached_pivot = this->get_pivot();
}

inline void VectorColumn::refresh(){
    return;
}




#endif // MPH_COLUMN_INCLUDED
