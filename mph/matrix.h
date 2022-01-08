#ifndef MPH_MATRIX_INCLUDED
#define MPH_MATRIX_INCLUDED

#include "signatureColumn.h"

struct Matrix : std::vector<SignatureColumn>{
    /** Encodes a matrix as a list of columns represented as lazy heaps. **/
    
    std::string to_string(){
        std::string output;
        hash_map<size_t, hash_set<size_t>> h_matrix;
        size_t max_row = 0;
        for(size_t column_nr = 0; column_nr < this->size(); column_nr++){
            for(auto& entry : this->at(column_nr)){
                size_t row_nr = entry;
                if(h_matrix.find(row_nr) == h_matrix.end()){
                    h_matrix[row_nr] = hash_set<size_t>();
                }
                if(h_matrix[row_nr].find(column_nr) != h_matrix[row_nr].end()){
                    h_matrix[row_nr].erase(column_nr);
                }else{
                    h_matrix[row_nr].insert(column_nr);
                }
                if(row_nr > max_row){
                    max_row = row_nr;
                }
            }
        }
        for(size_t row_nr=0; row_nr <= max_row; row_nr++){
            if(h_matrix.find(row_nr) == h_matrix.end()){
                for(size_t i=0; i<this->size(); i++){
                    output += "0 ";
                }
                output += "\n";
            }else{
                for(size_t i=0; i<this->size(); i++){
                    if(h_matrix[row_nr].find(i) == h_matrix[row_nr].end()){
                        output += "0 ";
                    }else{
                        output += "1 ";
                    }
                }
                output += "\n";
            }
        }
        return output;
    }
    
    void print(){
        std::cout << to_string();
    }
    
};

template <typename Column> void reduce_columns_vect2(std::vector<Column>& columns, index_t index, std::vector<index_t>& pivot_map) {
    /**  The function reduces a subset of the columns 'columns' described by the indices 'indices'. The function returns a pair containing a list of reduced columns and their syzygy-columns describing the operations performed on the columns. Only columns that were modified by the function will be returned in this list.
     
     Arguments:
     columns {std::vector<SignatureColumn>} -- a list of columns.
     indices {std::vector<index_t>} -- a list of indices that specify a submatrix from the list 'columns' that should be reduced.
     
     Returns:
     std::vector<SignatureColumn> -- a list of columns that were modified by the function in reduced form.
     */
    index_t pivot = columns[index].get_pivot().get_index();
    if(pivot_map[pivot] > -1 && pivot_map[pivot] < index){
        Column& working_column = columns[index];
        while(pivot != -1 && pivot_map[pivot] > -1 && pivot_map[pivot] < index){
            working_column.plus(columns[pivot_map[pivot]]);
            pivot = working_column.get_pivot().get_index();
        }
        if(pivot != -1){
            pivot_map[pivot] = index;
            working_column.refresh();
            working_column.syzygy.refresh();
        }
    } else{
        if(pivot != -1){
            pivot_map[pivot] = index;
        }
    }
    return;
}

template <typename Column> void reduce_columns_vect(std::vector<Column>& columns, std::vector<index_t>& indices, std::vector<Column>& modified_columns, std::vector<index_t>& pivot_map) {
    /**  The function reduces a subset of the columns 'columns' described by the indices 'indices'. The function returns a pair containing a list of reduced columns and their syzygy-columns describing the operations performed on the columns. Only columns that were modified by the function will be returned in this list.
     
     Arguments:
     columns {std::vector<SignatureColumn>} -- a list of columns.
     indices {std::vector<index_t>} -- a list of indices that specify a submatrix from the list 'columns' that should be reduced.
     
     Returns:
     std::vector<SignatureColumn> -- a list of columns that were modified by the function in reduced form.
     */
    
    
    for(size_t i=0; i<indices.size(); i++){
        if(pivot_map[columns[indices[i]].get_pivot().get_index()] > -1){
            index_t pivot = columns[indices[i]].get_pivot().get_index();
            Column working_column(columns[indices[i]]);
            
            int j=0;
            while(pivot != -1 && pivot_map[pivot] > -1){
                working_column.plus(columns[pivot_map[pivot]]);
                pivot = working_column.get_pivot().get_index();
                j++;
            }
            
            if(!working_column.empty()){
                pivot_map[working_column.get_pivot().get_index()] = indices[i];
            }
            
            modified_columns.push_back(working_column);
        } else{
            if(!columns[indices[i]].empty()){
                pivot_map[columns[indices[i]].get_pivot().get_index()] = indices[i];
            }
        }
    }
    return;
}

template <typename Column> void reduce_columns_temp(std::vector<Column>& columns, std::vector<index_t>& indices, std::vector<Column>& modified_columns, hash_map<index_t, index_t>& pivot_map) {
    /**  The function reduces a subset of the columns 'columns' described by the indices 'indices'. The function returns a pair containing a list of reduced columns and their syzygy-columns describing the operations performed on the columns. Only columns that were modified by the function will be returned in this list.
     
     Arguments:
     columns {std::vector<SignatureColumn>} -- a list of columns.
     indices {std::vector<index_t>} -- a list of indices that specify a submatrix from the list 'columns' that should be reduced.
     
     Returns:
     std::vector<SignatureColumn> -- a list of columns that were modified by the function in reduced form.
     */
    
    
    for(size_t i=0; i<indices.size(); i++){
        if(pivot_map.find(columns[indices[i]].get_pivot().get_index()) != pivot_map.end()){
            index_t pivot = columns[indices[i]].get_pivot().get_index();
            Column working_column(columns[indices[i]]);
            
            while(!working_column.empty() && pivot_map.find(pivot) != pivot_map.end()){
                working_column.plus(columns[pivot_map[pivot]]);
                pivot = working_column.get_pivot().get_index();
            }
            
            if(!working_column.empty()){
                pivot_map[working_column.get_pivot().get_index()] = indices[i];
            }
            
            modified_columns.push_back(working_column);
        } else{
            if(!columns[indices[i]].empty()){
                pivot_map[columns[indices[i]].get_pivot().get_index()] = indices[i];
            }
        }
    }
    return;
}

void reduce_columns(std::vector<SignatureColumn>& columns, std::vector<index_t>& indices, std::vector<SignatureColumn>& modified_columns, bool control_size=false) {
    /**  The function reduces a subset of the columns 'columns' described by the indices 'indices'. The function returns a pair containing a list of reduced columns and their syzygy-columns describing the operations performed on the columns. Only columns that were modified by the function will be returned in this list.
     
     Arguments:
     columns {std::vector<SignatureColumn>} -- a list of columns.
     indices {std::vector<index_t>} -- a list of indices that specify a submatrix from the list 'columns' that should be reduced.
     
     Returns:
     std::vector<SignatureColumn> -- a list of columns that were modified by the function in reduced form.
     */
    
    hash_map<index_t, index_t> pivot_map;
    
    
    for(size_t i=0; i<indices.size(); i++){
        if(pivot_map.find(columns[indices[i]].get_pivot().get_index()) != pivot_map.end()){
            index_t pivot = columns[indices[i]].get_pivot().get_index();
            SignatureColumn working_column(columns[indices[i]].grade, columns[indices[i]].signature_index);
            working_column.plus(columns[indices[i]]);
            size_t input_size=working_column.size();

            while(!working_column.empty() && pivot_map.find(pivot) != pivot_map.end()){
                working_column.plus(columns[pivot_map[pivot]]);
                pivot = working_column.get_pivot().get_index();
                if(control_size && working_column.size() > 2*input_size){
                    working_column.refresh();
                    working_column.syzygy.refresh();
                    input_size = working_column.size();
                }
            }
            
            if(!working_column.empty()){
                pivot_map[working_column.get_pivot().get_index()] = indices[i];
            }
            
            modified_columns.push_back(working_column);
        } else{
            if(!columns[indices[i]].empty()){
                pivot_map[columns[indices[i]].get_pivot().get_index()] = indices[i];
            }
        }
    }
    return;
}




class CustomHashMap{
public:
    std::vector<std::vector<size_t>> keys;
    std::vector<std::vector<size_t>> values;
    
    CustomHashMap(size_t size){
        keys.reserve(2*size+1);
        values.reserve(2*size+1);
        for(size_t i=0; i<2*size+1; i++){
            keys.push_back(std::vector<size_t>());
            values.push_back(std::vector<size_t>());
        }
    }
    
    bool containsKey(size_t key){
        if(keys[key%keys.size()].size() == 0){
            return false;
        }else{
            for(auto entry : keys[key%keys.size()]){
                if(entry == key){
                    return true;
                }
            }
            return false;
        }
    }
    
    void insert(size_t key, size_t val){
        keys[key%keys.size()].push_back(key);
        values[key%keys.size()].push_back(val);
    }
    
    size_t get(size_t key){
        for(size_t i=0; i<keys[key%keys.size()].size(); i++){
            if(keys[key%keys.size()][i] == key){
                return values[key%keys.size()][i];
            }
        }
        throw "Could not find the values corresponding to key: "+std::to_string(key);
        return -1;
    }
    
};

void reduce_columns_custom_hash(Matrix& columns, std::vector<index_t>& indices, Matrix& modified_columns) {
    /**  The function reduces a subset of the columns 'columns' described by the indices 'indices'. The function returns a pair containing a list of reduced columns and their syzygy-columns describing the operations performed on the columns. Only columns that were modified by the function will be returned in this list.
     
     Arguments:
     columns {std::vector<SignatureColumn>} -- a list of columns.
     indices {std::vector<index_t>} -- a list of indices that specify a submatrix from the list 'columns' that should be reduced.
     
     Returns:
     std::vector<SignatureColumn> -- a list of columns that were modified by the function in reduced form.
     */
    
    CustomHashMap map(indices.size());
    
    
    for(size_t i=0; i<indices.size(); i++){
        if(map.containsKey(columns[indices[i]].get_pivot().get_index())){
            index_t pivot = columns[indices[i]].get_pivot().get_index();
            SignatureColumn working_column(columns[indices[i]].grade, columns[indices[i]].signature_index);
            working_column.plus(columns[indices[i]]);
            
            
            while(!working_column.empty() && map.containsKey(pivot)){
                working_column.plus(columns[map.get(pivot)]);
                pivot = working_column.get_pivot().get_index();
            }
            if(!working_column.empty()){
                map.insert(working_column.get_pivot().get_index(), indices[i]);
            }
            
            modified_columns.push_back(working_column);
        } else{
            if(!columns[indices[i]].empty()){
                map.insert(columns[indices[i]].get_pivot().get_index(),  indices[i]);
            }
        }
    }
    return;
}


void reduce_columns_par(Matrix& columns, Matrix& additional_columns, std::vector<index_t>& indices, std::vector<SignatureColumn*>& modified_columns) {
    /**  The function reduces a subset of the columns 'columns' described by the indices 'indices'. The function returns a pair containing a list of reduced columns and their syzygy-columns describing the operations performed on the columns. Only columns that were modified by the function will be returned in this list.
     
     Arguments:
     columns {std::vector<SignatureColumn>} -- a list of columns.
     indices {std::vector<index_t>} -- a list of indices that specify a submatrix from the list 'columns' that should be reduced.
     
     Returns:
     std::vector<SignatureColumn> -- a list of columns that were modified by the function in reduced form.
     */
    
    hash_map<index_t, index_t> pivot_map;
    
    
    
    for(size_t i=0; i<indices.size(); i++){
        bool base_list = indices[i] < columns.size();
        if(pivot_map.find((base_list ? columns[indices[i]] : additional_columns[indices[i]-columns.size()]).get_cached_pivot().get_index()) != pivot_map.end()){
            SignatureColumn working_column((base_list ? columns[indices[i]] : additional_columns[indices[i]-columns.size()]).grade, (base_list ? columns[indices[i]] : additional_columns[indices[i]-columns.size()]).signature_index);
            working_column.plus((base_list ? columns[indices[i]] : additional_columns[indices[i]-columns.size()]));
            index_t pivot = working_column.get_pivot().get_index();
            
            while(!working_column.empty() && (pivot_map.find(pivot) != pivot_map.end())){
                working_column.plus(pivot_map[pivot] < columns.size() ? columns[pivot_map[pivot]] : additional_columns[pivot_map[pivot]-columns.size()]);
                pivot = working_column.get_pivot().get_index();
            }
            if(!working_column.empty()){
                pivot_map[working_column.get_pivot().get_index()] = indices[i];
            }
            modified_columns.push_back(new SignatureColumn(working_column));
        } else{
            if((base_list ? columns[indices[i]] : additional_columns[indices[i]-columns.size()]).get_cached_pivot().get_index() != -1){
                pivot_map[(base_list ? columns[indices[i]] : additional_columns[indices[i]-columns.size()]).get_cached_pivot().get_index()] = indices[i];
            }
        }
    }
    return;
}






#endif // MPH_MATRIX_INCLUDED
