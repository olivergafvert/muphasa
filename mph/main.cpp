//
//  main.cpp
//  mph
//
//  Created by Oliver on 2020-04-09.
//  Copyright © 2020 Oliver. All rights reserved.
//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
#include <queue>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <set>
#include <random>
#include <time.h>
#include <boost/progress.hpp>

#include "utils.h"
#include "grade.h"
#include "column.h"
#include "signatureColumn.h"
#include "matrix.h"
#include "IO.h"
#include "examples.h"


double res_memory=0;
double virt_memory=0;

/* Groebner bases */

std::pair<Matrix, Matrix> computeGroebnerBases(std::vector<SignatureColumn>& columns){
    /*
     The main function computing a Groebner basis for the image and kernel of the map described by the list of columns 'columns'.
     
     Arguments:
     columns {std::vector<SignatureColumn>} -- columns describing the matrix of a map between two free multigraded momdules.
     
     Returns:
     std::vector<SignatureColumn> -- a list of vectors decribing a minimal Groebner basis for the image of the map.
     std::vector<SignatureColumn> -- a list of vectors describing a minimal Groebner basis for the kernel of the map.
     */
    
    std::cout << "Starting to compute Groebner bases..." << std::endl;
    
    /* Sort columns colexicographically */
    sort(columns.begin(), columns.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });
    // The sorted columns should agree with the columns sorted by index of signature
    hash_map<size_t, size_t> index_map_high;
    for(size_t i=0; i<columns.size(); i++){
        columns[i].signature_index = i;
        index_map_high[columns[i].grade[columns[i].grade.size()-1]] = i;
    }
    
    /* Compute index set iterator */
    std::vector<std::vector<index_t>> grade_base_set = get_grade_base_set<SignatureColumn>(columns);
    int n_total_grades = 1;
    for(auto& grade_list : grade_base_set){
        n_total_grades *= grade_list.size();
    }
    Iterator_lex grade_iterator = Iterator_lex(grade_base_set);
    
    /* Vectors to store the columns of the GBs */
    std::vector<SignatureColumn> gb_columns;
    std::vector<SignatureColumn> syzygies;
    
    /* Index maps to keep track of signatures */
    std::vector<std::vector<signature_t>> GB;
    std::vector<std::vector<signature_t>> Syz;
    
    for(size_t i=0; i<columns.size(); i++){
        GB.push_back(std::vector<signature_t>());
        Syz.push_back(std::vector<signature_t>());
    }
    
    /* Main algorithm that iterates through the index set */
    int column_index;
    index_t max_pivot=0;
    for(auto& column : columns){
        if(max_pivot < column.get_pivot().get_index()){
            max_pivot = column.get_pivot().get_index();
        }
    }
    index_t pivot;
    
    int iter_index = 0;
    
    while(grade_iterator.has_next()){
        grade_t& v = grade_iterator.next();
        iter_index++;
        std::vector<index_t> pivot_map(max_pivot+1, -1);
        /* Initialize Macaulay matrix */
        size_t& index_bound = index_map_high[v[v.size()-1]];
        for( size_t i=0; i<=index_bound; i++ ){
            column_index = -1;
            bool is_new=false;
            if(GB[i].size() > 0){
                bool in_syz = false;
                if(Syz[i].size()>0){
                    for(size_t j=Syz[i].size()-1; j < Syz[i].size(); j--){
                        if((Syz[i][j].get_grade()).leq_poset(v)){
                            in_syz = true;
                            break;
                        }
                    }
                }
                if(!in_syz){
                    for(size_t j=GB[i].size()-1; j < GB[i].size(); j--){
                        if(GB[i][j].get_grade().leq_poset(v)){
                            column_index = (int)GB[i][j].get_index();
                            break;
                        }
                    }
                }
            } else{
                if(columns[i].grade == v){
                    gb_columns.push_back(columns[i]);
                    column_index = (int)gb_columns.size()-1;
                    is_new=true;
                }
            }
            if(column_index > -1){
                pivot = gb_columns[column_index].get_pivot().get_index();
                if(pivot_map[pivot] > -1){// && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                    SignatureColumn working_column = gb_columns[column_index];
                    while(pivot != -1 && pivot_map[pivot] > -1){// && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                        working_column.plus(gb_columns[pivot_map[pivot]]);
                        pivot = working_column.get_pivot().get_index();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        gb_columns.push_back(working_column);
                        GB[i].push_back(signature_t(working_column.grade, gb_columns.size()-1));
                        pivot_map[pivot] = gb_columns.size()-1;
                        gb_columns[gb_columns.size()-1].last_updated = iter_index;
                    }else{
                        working_column.syzygy.refresh();
                        syzygies.push_back(SignatureColumn(working_column.grade, working_column.signature_index, working_column.syzygy));
                        Syz[i].push_back(signature_t(syzygies.back().grade, syzygies.size()-1));
                    }
                } else{
                    if(pivot != -1){
                        pivot_map[pivot] = column_index;
                        gb_columns[column_index].last_updated = iter_index;
                        if(is_new){
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, column_index));
                        }
                    }
                }
            }
        }
    }
    
    Matrix syzygies_output;
    syzygies_output.reserve(syzygies.size());
    for(size_t i=0; i<Syz.size(); i++){
        for(size_t j=0; j<Syz[i].size(); j++){
            syzygies_output.push_back(SignatureColumn(Syz[i][j].get_grade(), syzygies_output.size(), syzygies[Syz[i][j].get_index()]));
        }
    }
    Matrix gb_columns_output;
    gb_columns_output.reserve(gb_columns.size());
    for(size_t i=0; i<GB.size(); i++){
        for(size_t j=0; j<GB[i].size(); j++){
            gb_columns_output.push_back(gb_columns[GB[i][j].get_index()]);
        }
    }
    std::cout << "Finished computing Groebner bases." << std::endl;
    return std::pair<Matrix, Matrix>(gb_columns_output, syzygies_output);
}

std::pair<Matrix, Matrix> computeGroebnerBases_gradeopt_min(std::vector<SignatureColumn>& columns){
    /*
     The main function computing a Groebner basis for the image and kernel of the map described by the list of columns 'columns'.
     
     Arguments:
     columns {std::vector<SignatureColumn>} -- columns describing the matrix of a map between two free multigraded momdules.
     
     Returns:
     std::vector<SignatureColumn> -- a list of vectors decribing a minimal Groebner basis for the image of the map.
     std::vector<SignatureColumn> -- a list of vectors describing a minimal Groebner basis for the kernel of the map.
     */
    
    std::cout << "Starting to compute Groebner bases..." << std::endl;
    
    /* Sort columns colexicographically */
    sort(columns.begin(), columns.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });
    // The sorted columns should agree with the columns sorted by index of signature
    hash_map<size_t, size_t> index_map_high;
    for(size_t i=0; i<columns.size(); i++){
        columns[i].signature_index = i;
        index_map_high[columns[i].grade[columns[i].grade.size()-1]] = i;
    }
    
    /* Compute index set iterator */
    std::priority_queue<grade_t, std::vector<grade_t>, std::greater<grade_t>> grades;
    std::vector<std::vector<grade_t>> grade_lists;
    std::unordered_set<grade_t, GradeHasher> visited_grades;
    
    /* Vectors to store the columns of the GBs */
    std::vector<SignatureColumn> gb_columns;
    std::vector<VectorColumn> syzygies;
    gb_columns.reserve(2*columns.size());
    syzygies.reserve(columns.size());
    
    /* Index maps to keep track of signatures */
    std::vector<std::vector<signature_t>> GB;
    std::vector<std::vector<signature_t>> Syz;
    GB.reserve(columns.size());
    Syz.reserve(columns.size());
    
    for(size_t i=0; i<columns.size(); i++){
        GB.push_back(std::vector<signature_t>());
        Syz.push_back(std::vector<signature_t>());
    }
    
    /* Main algorithm that iterates through the index set */
    int column_index;
    index_t max_pivot=0;
    for(auto& column : columns){
        if(max_pivot < column.get_pivot_index()){
            max_pivot = column.get_pivot_index();
        }
    }
    for(size_t i=0; i<=max_pivot; i++){
        grade_lists.push_back(std::vector<grade_t>());
    }
    for(auto& column : columns){
        if(visited_grades.find(column.grade) == visited_grades.end()){
            grades.push(column.grade);
            visited_grades.insert(column.grade);
        }
    }
    
    std::vector<size_t> grade_hashes;
    grade_hashes.reserve(columns.size());
    GradeHasher grade_hasher;
    for(auto& c : columns){
        grade_hashes.push_back(grade_hasher(c.grade));
    }
    std::vector<index_t> pivot_map(max_pivot+1, -1);
    index_t pivot;
    
    int iter_index = 0;
    
    while(!grades.empty()){
        grade_t v = grades.top();
        grades.pop();
        while(v == grades.top()){
            grades.pop();
        }
        
        size_t grade_hash = grade_hasher(v);
        iter_index++;
        
        /* Initialize Macaulay matrix */
        size_t& index_bound = index_map_high[v[v.size()-1]];
        for( size_t i=0; i<=index_bound; i++ ){
            column_index = -1;
            bool is_new=false;
            if(GB[i].size() > 0){
                bool in_syz = false;
                if(Syz[i].size()>0){
                    for(size_t j=Syz[i].size()-1; j < Syz[i].size(); j--){
                        if((Syz[i][j].get_grade()).leq_poset(v)){
                            in_syz = true;
                            break;
                        }
                    }
                }
                if(!in_syz){
                    for(size_t j=GB[i].size()-1; j < GB[i].size(); j--){
                        if(GB[i][j].get_grade().leq_poset(v)){
                            column_index = (int)GB[i][j].get_index();
                            break;
                        }
                    }
                }
            } else{
                if(grade_hash == grade_hashes[i] && columns[i].grade == v){
                    gb_columns.push_back(columns[i]);
                    column_index = (int)gb_columns.size()-1;
                    is_new=true;
                }
            }
            if(column_index > -1){
                pivot = gb_columns[column_index].get_pivot_index();
                if(pivot_map[pivot] > -1 && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                    SignatureColumn working_column = gb_columns[column_index];
                    while(pivot != -1 && pivot_map[pivot] > -1 && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                        working_column.plus(gb_columns[pivot_map[pivot]]);
                        pivot = working_column.get_pivot_index();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        gb_columns.push_back(working_column);
                        GB[i].push_back(signature_t(working_column.grade, gb_columns.size()-1));
                        pivot_map[pivot] = gb_columns.size()-1;
                        gb_columns[gb_columns.size()-1].last_updated = iter_index;
                        if(working_column.grade == v){
                            if(grade_lists[pivot].size() == 0 || working_column.grade != grade_lists[pivot].back()){
                                std::vector<grade_t> minimal_elements_tmp;
                                for(size_t j=0; j<grade_lists[pivot].size(); j++){
                                    grade_t m_ji = grade_lists[pivot][j].m_ji(working_column.grade);
                                    bool is_minimal = true;
                                    for(auto& el : minimal_elements_tmp){
                                        if(el.leq_poset(m_ji)){
                                            is_minimal = false;
                                            break;
                                        }
                                    }
                                    if(is_minimal){
                                        minimal_elements_tmp.push_back(m_ji);
                                        grade_t g = working_column.grade.join(grade_lists[pivot][j]);
                                        if(visited_grades.find(g) == visited_grades.end()){
                                            grades.push(g);
                                            visited_grades.insert(g);
                                        }
                                    }
                                }
                                grade_lists[pivot].push_back(working_column.grade);
                            }
                        }
                    }else{
                        working_column.syzygy.refresh();
                        syzygies.push_back(working_column.syzygy);
                        Syz[i].push_back(signature_t(working_column.grade, syzygies.size()-1));
                    }
                } else{
                    if(pivot != -1){
                        pivot_map[pivot] = column_index;
                        gb_columns[column_index].last_updated = iter_index;
                        if(is_new){
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, column_index));
                            if(grade_lists[pivot].size() == 0 || gb_columns[column_index].grade != grade_lists[pivot].back()){
                                std::vector<grade_t> minimal_elements_tmp;
                                for(size_t j=0; j<grade_lists[pivot].size(); j++){
                                    grade_t m_ji = grade_lists[pivot][j].m_ji(gb_columns[column_index].grade);
                                    bool is_minimal = true;
                                    for(auto& el : minimal_elements_tmp){
                                        if(el.leq_poset(m_ji)){
                                            is_minimal = false;
                                            break;
                                        }
                                    }
                                    if(is_minimal){
                                        minimal_elements_tmp.push_back(m_ji);
                                        grade_t g = gb_columns[column_index].grade.join(grade_lists[pivot][j]);
                                        if(visited_grades.find(g) == visited_grades.end()){
                                            grades.push(g);
                                            visited_grades.insert(g);
                                        }
                                    }
                                }
                                grade_lists[pivot].push_back(gb_columns[column_index].grade);
                            }
                        }
                    }
                }
            }
        }
    }
    
    Matrix syzygies_output;
    syzygies_output.reserve(syzygies.size());
    for(size_t i=0; i<Syz.size(); i++){
        for(size_t j=0; j<Syz[i].size(); j++){
            syzygies_output.push_back(SignatureColumn(Syz[i][j].get_grade(), syzygies_output.size(), syzygies[Syz[i][j].get_index()]));
        }
    }
    Matrix gb_columns_output;
    gb_columns_output.reserve(gb_columns.size());
    for(size_t i=0; i<GB.size(); i++){
        for(size_t j=0; j<GB[i].size(); j++){
            gb_columns_output.push_back(gb_columns[GB[i][j].get_index()]);
        }
    }
    std::cout << "Finished computing Groebner bases." << std::endl;
    return std::pair<Matrix, Matrix>(gb_columns_output, syzygies_output);
}

std::pair<Matrix, Matrix> computeGroebnerBases_gradeopt(std::vector<SignatureColumn>& columns){
    /*
     The main function computing a Groebner basis for the image and kernel of the map described by the list of columns 'columns'.
     
     Arguments:
     columns {std::vector<SignatureColumn>} -- columns describing the matrix of a map between two free multigraded momdules.
     
     Returns:
     std::vector<SignatureColumn> -- a list of vectors decribing a minimal Groebner basis for the image of the map.
     std::vector<SignatureColumn> -- a list of vectors describing a minimal Groebner basis for the kernel of the map.
     */
    
    std::cout << "Starting to compute Groebner bases..." << std::endl;
    
    /* Sort columns colexicographically */
    sort(columns.begin(), columns.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });
    // The sorted columns should agree with the columns sorted by index of signature
    hash_map<size_t, size_t> index_map_high;
    for(size_t i=0; i<columns.size(); i++){
        columns[i].signature_index = i;
        index_map_high[columns[i].grade[columns[i].grade.size()-1]] = i;
    }
    
    /* Compute index set iterator */
    std::priority_queue<grade_t, std::vector<grade_t>, std::greater<grade_t>> grades;
    std::vector<std::vector<grade_t>> grade_lists;
    std::unordered_set<grade_t, GradeHasher> visited_grades;
    
    /* Vectors to store the columns of the GBs */
    Matrix gb_columns;
    Matrix syzygies;
    gb_columns.reserve(columns.size());
    syzygies.reserve(columns.size());
    
    /* Index maps to keep track of signatures */
    std::vector<std::vector<signature_t>> GB;
    std::vector<std::vector<signature_t>> Syz;
    GB.reserve(columns.size());
    Syz.reserve(columns.size());
    
    for(size_t i=0; i<columns.size(); i++){
        GB.push_back(std::vector<signature_t>());
        Syz.push_back(std::vector<signature_t>());
    }
    
    /* Main algorithm that iterates through the index set */
    int column_index;
    index_t max_pivot=0;
    for(auto& column : columns){
        if(max_pivot < column.get_pivot_index()){
            max_pivot = column.get_pivot_index();
        }
    }
    for(size_t i=0; i<=max_pivot; i++){
        grade_lists.push_back(std::vector<grade_t>());
    }
    for(auto& column : columns){
        grades.push(column.grade);
    }
    
    std::vector<size_t> grade_hashes;
    grade_hashes.reserve(columns.size());
    GradeHasher grade_hasher;
    for(auto& c : columns){
        grade_hashes.push_back(grade_hasher(c.grade));
    }
    
    index_t pivot;
    std::vector<index_t> pivot_map(max_pivot+1, -1);
    int iter_index = 0;
    
    while(!grades.empty()){
        grade_t v = grades.top();
        grades.pop();
        while(v == grades.top()){
            grades.pop();
        }
        size_t grade_hash = grade_hasher(v);
        iter_index++;
        
        /* Initialize Macaulay matrix */
        size_t& index_bound = index_map_high[v[v.size()-1]];
        for( size_t i=0; i<=index_bound; i++ ){
            column_index = -1;
            bool is_new=false;
            if(GB[i].size() > 0){
                bool in_syz = false;
                if(Syz[i].size()>0){
                    for(size_t j=Syz[i].size()-1; j < Syz[i].size(); j--){
                        if((Syz[i][j].get_grade()).leq_poset(v)){
                            in_syz = true;
                            break;
                        }
                    }
                }
                if(!in_syz){
                    for(size_t j=GB[i].size()-1; j < GB[i].size(); j--){
                        if(GB[i][j].get_grade().leq_poset(v)){
                            column_index = (int)GB[i][j].get_index();
                            break;
                        }
                    }
                }
            } else{
                if(grade_hash == grade_hashes[i] && columns[i].grade == v){
                    gb_columns.push_back(columns[i]);
                    column_index = (int)gb_columns.size()-1;
                    is_new=true;
                }
            }
            if(column_index > -1){
                pivot = gb_columns[column_index].get_pivot_index();
                if(pivot_map[pivot] != -1 && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                    SignatureColumn working_column = gb_columns[column_index];
                    while(pivot != -1 && pivot_map[pivot] > -1 && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                        working_column.plus(gb_columns[pivot_map[pivot]]);
                        pivot = working_column.get_pivot_index();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        grade_t grade = working_column.grade;
                        if(is_new){
                            gb_columns[column_index] = working_column;
                            GB[i].push_back(signature_t(grade, column_index));
                            pivot_map[pivot] = column_index;
                            gb_columns[column_index].last_updated = iter_index;
                        }else{
                            /*if(GB[i].size()==1){
                                gb_columns.push_back(working_column);
                                GB[i].push_back(signature_t(grade, gb_columns.size()-1));
                            }else{
                                if(pivot_map[gb_columns[GB[i][1].get_index()].get_pivot_index()] == GB[i][1].get_index()){
                                    pivot_map[gb_columns[GB[i][1].get_index()].get_pivot_index()] = -1;
                                }
                                gb_columns[GB[i][1].get_index()].swap(working_column);
                                GB[i][1].first = grade;
                            }*/
                            gb_columns.push_back(working_column);
                            GB[i].push_back(signature_t(grade, gb_columns.size()-1));
                            pivot_map[pivot] = gb_columns.size()-1;
                            gb_columns[gb_columns.size()-1].last_updated = iter_index;
                        }
                        if(grade == v){
                            if(grade_lists[pivot].size() == 0 || grade != grade_lists[pivot].back()){
                                for(size_t ig=0; ig<grade_lists[pivot].size(); ig++){
                                    grade_t g = grade.join(grade_lists[pivot][ig]);
                                    if(grade != g){
                                        grades.push(g);
                                    }
                                }
                                grade_lists[pivot].push_back(grade);
                            }
                        }
                       
                    }else{
                        if(is_new){
                            gb_columns.pop_back();
                        }
                        working_column.syzygy.refresh();
                        syzygies.push_back(SignatureColumn(working_column.get_grade(), syzygies.size(), working_column.syzygy));
                        Syz[i].push_back(signature_t(working_column.grade, syzygies.size()-1));
                    }
                } else{
                    if(pivot != -1){
                        pivot_map[pivot] = column_index;
                        gb_columns[column_index].last_updated = iter_index;
                        if(is_new){
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, column_index));
                            if(grade_lists[pivot].size() == 0 || gb_columns[column_index].grade != grade_lists[pivot].back()){
                                for(size_t ig=0; ig<grade_lists[pivot].size(); ig++){
                                    grade_t g = gb_columns[column_index].grade.join(grade_lists[pivot][ig]);
                                    if(gb_columns[column_index].grade != g){
                                        grades.push(g);
                                    }
                                }
                                grade_lists[pivot].push_back(gb_columns[column_index].grade);
                            }
                        }
                    }
                }
            }
        }
    }
    
    get_mem_usage(virt_memory, res_memory);
    
    std::cout << "Finished computing Groebner bases." << std::endl;
    return std::pair<Matrix, Matrix>(gb_columns, syzygies);
}

Matrix computekernel_gradeopt(std::vector<SignatureColumn>& columns){
    /*
     The main function computing a Groebner basis for the image and kernel of the map described by the list of columns 'columns'.
     
     Arguments:
     columns {std::vector<SignatureColumn>} -- columns describing the matrix of a map between two free multigraded momdules.
     
     Returns:
     std::vector<SignatureColumn> -- a list of vectors decribing a minimal Groebner basis for the image of the map.
     std::vector<SignatureColumn> -- a list of vectors describing a minimal Groebner basis for the kernel of the map.
     */
    
    std::cout << "Starting to compute Groebner bases..." << std::endl;
    
    /* Sort columns colexicographically */
    sort(columns.begin(), columns.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });
    // The sorted columns should agree with the columns sorted by index of signature
    hash_map<size_t, size_t> index_map_high;
    for(size_t i=0; i<columns.size(); i++){
        columns[i].signature_index = i;
        index_map_high[columns[i].grade[columns[i].grade.size()-1]] = i;
    }
    
    /* Compute index set iterator */
    std::priority_queue<grade_t, std::vector<grade_t>, std::greater<grade_t>> grades;
    std::vector<std::vector<grade_t>> grade_lists;
    
    /* Vectors to store the columns of the GBs */
    std::vector<SignatureColumn> gb_columns;
    std::vector<SyzColumn> syzygies;
    gb_columns.reserve(2*columns.size());
    syzygies.reserve(columns.size());
    
    /* Index maps to keep track of signatures */
    std::vector<std::vector<signature_t>> GB;
    std::vector<std::vector<signature_t>> Syz;
    GB.reserve(columns.size());
    Syz.reserve(columns.size());
    
    for(size_t i=0; i<columns.size(); i++){
        GB.push_back(std::vector<signature_t>());
        Syz.push_back(std::vector<signature_t>());
    }
    
    /* Main algorithm that iterates through the index set */
    int column_index;
    index_t max_pivot=0;
    for(auto& column : columns){
        if(max_pivot < column.get_pivot_index()){
            max_pivot = column.get_pivot_index();
        }
    }
    for(size_t i=0; i<=max_pivot; i++){
        grade_lists.push_back(std::vector<grade_t>());
    }
    for(auto& column : columns){
        grades.push(column.grade);
    }
    
    index_t pivot;
    
    int iter_index = 0;
    
    while(!grades.empty()){
        grade_t v = grades.top();
        grades.pop();
        while(v == grades.top()){
            grades.pop();
        }
        iter_index++;
        std::vector<index_t> pivot_map(max_pivot+1, -1);
        /* Initialize Macaulay matrix */
        size_t& index_bound = index_map_high[v[v.size()-1]];
        for( size_t i=0; i<=index_bound; i++ ){
            column_index = -1;
            bool is_new=false;
            if(GB[i].size() > 0){
                bool in_syz = false;
                if(Syz[i].size()>0){
                    for(size_t j=Syz[i].size()-1; j < Syz[i].size(); j--){
                        if((Syz[i][j].get_grade()).leq_poset(v)){
                            in_syz = true;
                            break;
                        }
                    }
                }
                if(!in_syz){
                    for(size_t j=GB[i].size()-1; j < GB[i].size(); j--){
                        if(GB[i][j].get_grade().leq_poset(v)){
                            column_index = (int)GB[i][j].get_index();
                            break;
                        }
                    }
                }
            } else{
                if(columns[i].grade == v){
                    gb_columns.push_back(columns[i]);
                    column_index = (int)gb_columns.size()-1;
                    is_new=true;
                }
            }
            if(column_index > -1){
                pivot = gb_columns[column_index].get_pivot_index();
                if(pivot_map[pivot] != -1){// && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                    SignatureColumn working_column = gb_columns[column_index];
                    while(pivot != -1 && pivot_map[pivot] > -1){// && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                        working_column.plus(gb_columns[pivot_map[pivot]]);
                        pivot = working_column.get_pivot_index();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        if(is_new){
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, column_index));
                            gb_columns.push_back(working_column);
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, gb_columns.size()-1));
                            pivot_map[pivot] = gb_columns.size()-1;
                            gb_columns[gb_columns.size()-1].last_updated = iter_index;
                        }else{
                            GB[i][1].first = working_column.grade;
                            gb_columns[GB[i][1].get_index()].swap(working_column);
                            pivot_map[pivot] = GB[i][1].get_index();
                            gb_columns[GB[i][1].get_index()].last_updated = iter_index;
                        }
                        grade_t& grade = GB[i][1].get_grade();
                        if(grade == v){
                            if(grade_lists[pivot].size() == 0 || grade != grade_lists[pivot].back()){
                                for(size_t ig=0; ig<grade_lists[pivot].size(); ig++){
                                    if(grade.join(grade_lists[pivot][ig]) != grade){
                                        grade_t g = grade_lists[pivot][ig];
                                        grades.push(grade.join(grade_lists[pivot][ig]));
                                    }
                                }
                                grade_lists[pivot].push_back(grade);
                            }
                        }
                       
                    }else{
                        working_column.syzygy.refresh();
                        syzygies.push_back(working_column.syzygy);
                        Syz[i].push_back(signature_t(working_column.grade, syzygies.size()-1));
                    }
                } else{
                    if(pivot != -1){
                        pivot_map[pivot] = column_index;
                        gb_columns[column_index].last_updated = iter_index;
                        if(is_new){
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, column_index));
                            gb_columns.push_back(gb_columns[column_index]);
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, gb_columns.size()-1));
                            pivot_map[pivot] = gb_columns.size()-1;
                            gb_columns[gb_columns.size()-1].last_updated = iter_index;
                            if(grade_lists[pivot].size() == 0 || gb_columns[column_index].grade != grade_lists[pivot].back()){
                                for(size_t ig=0; ig<grade_lists[pivot].size(); ig++){
                                    if(gb_columns[column_index].grade.join(grade_lists[pivot][ig]) != gb_columns[column_index].grade){
                                        grades.push(gb_columns[column_index].grade.join(grade_lists[pivot][ig]));
                                    }
                                }
                                grade_lists[pivot].push_back(gb_columns[column_index].grade);
                            }
                        }
                    }
                }
            }
        }
    }
    
    Matrix syzygies_output;
    syzygies_output.reserve(syzygies.size());
    for(size_t i=0; i<Syz.size(); i++){
        for(size_t j=0; j<Syz[i].size(); j++){
            syzygies_output.push_back(SignatureColumn(Syz[i][j].get_grade(), syzygies_output.size(), syzygies[Syz[i][j].get_index()]));
        }
    }
    std::cout << "Finished computing Groebner bases." << std::endl;
    return syzygies_output;
}

Matrix ImageGB(std::vector<SignatureColumn>& columns){
    /*
     The main function computing a Groebner basis for the image and kernel of the map described by the list of columns 'columns'.
     
     Arguments:
     columns {std::vector<SignatureColumn>} -- columns describing the matrix of a map between two free multigraded momdules.
     
     Returns:
     std::vector<SignatureColumn> -- a list of vectors decribing a minimal Groebner basis for the image of the map.
     std::vector<SignatureColumn> -- a list of vectors describing a minimal Groebner basis for the kernel of the map.
     */
    
    std::cout << "Starting to compute Groebner bases..." << std::endl;
    
    /* Sort columns colexicographically */
    sort(columns.begin(), columns.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });
    // The sorted columns should agree with the columns sorted by index of signature
    hash_map<size_t, size_t> index_map_high;
    for(size_t i=0; i<columns.size(); i++){
        columns[i].signature_index = i;
        index_map_high[columns[i].grade[columns[i].grade.size()-1]] = i;
    }
    
    /* Compute index set iterator */
    std::priority_queue<grade_t, std::vector<grade_t>, std::greater<grade_t>> grades;
    std::vector<std::vector<grade_t>> grade_lists;
    std::unordered_set<grade_t, GradeHasher> visited_grades;
    
    /* Vectors to store the columns of the GBs */
    Matrix gb_columns;
    gb_columns.reserve(columns.size());
    
    /* Index maps to keep track of signatures */
    std::vector<std::vector<signature_t>> GB;
    std::vector<std::vector<signature_t>> Syz;
    GB.reserve(columns.size());
    Syz.reserve(columns.size());
    
    for(size_t i=0; i<columns.size(); i++){
        GB.push_back(std::vector<signature_t>());
        Syz.push_back(std::vector<signature_t>());
    }
    
    /* Main algorithm that iterates through the index set */
    int column_index;
    index_t max_pivot=0;
    for(auto& column : columns){
        if(max_pivot < column.get_pivot_index()){
            max_pivot = column.get_pivot_index();
        }
    }
    for(size_t i=0; i<=max_pivot; i++){
        grade_lists.push_back(std::vector<grade_t>());
    }
    for(auto& column : columns){
        if(visited_grades.find(column.grade) == visited_grades.end()){
            grades.push(column.grade);
            visited_grades.insert(column.grade);
        }
    }
    
    std::vector<size_t> grade_hashes;
    grade_hashes.reserve(columns.size());
    GradeHasher grade_hasher;
    for(auto& c : columns){
        grade_hashes.push_back(grade_hasher(c.grade));
    }
    
    index_t pivot;
    std::vector<index_t> pivot_map(max_pivot+1, -1);
    int iter_index = 0;
    
    while(!grades.empty()){
        grade_t v = grades.top();
        grades.pop();
        /*while(v == grades.top()){
            grades.pop();
        }*/
        size_t grade_hash = grade_hasher(v);
        iter_index++;
        
        /* Initialize Macaulay matrix */
        size_t& index_bound = index_map_high[v[v.size()-1]];
        for( size_t i=0; i<=index_bound; i++ ){
            column_index = -1;
            bool is_new=false;
            if(GB[i].size() > 0){
                bool in_syz = false;
                if(Syz[i].size()>0){
                    for(size_t j=Syz[i].size()-1; j < Syz[i].size(); j--){
                        if((Syz[i][j].get_grade()).leq_poset(v)){
                            in_syz = true;
                            break;
                        }
                    }
                }
                if(!in_syz){
                    for(size_t j=GB[i].size()-1; j < GB[i].size(); j--){
                        if(GB[i][j].get_grade().leq_poset(v)){
                            column_index = (int)GB[i][j].get_index();
                            break;
                        }
                    }
                }
            } else{
                if(grade_hash == grade_hashes[i] && columns[i].grade == v){
                    gb_columns.push_back(columns[i]);
                    column_index = (int)gb_columns.size()-1;
                    is_new=true;
                }
            }
            if(column_index > -1){
                pivot = gb_columns[column_index].get_pivot_index();
                if(pivot_map[pivot] != -1 && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                    SignatureColumn working_column = gb_columns[column_index];
                    while(pivot != -1 && pivot_map[pivot] > -1 && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                        working_column.plus(gb_columns[pivot_map[pivot]]);
                        pivot = working_column.get_pivot_index();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        if(is_new){
                            gb_columns.pop_back();
                        }
                        gb_columns.push_back(working_column);
                        GB[i].push_back(signature_t(working_column.grade, gb_columns.size()-1));
                        pivot_map[pivot] = gb_columns.size()-1;
                        gb_columns[gb_columns.size()-1].last_updated = iter_index;
                        if(working_column.grade == v){
                            if(grade_lists[pivot].size() == 0 || working_column.grade != grade_lists[pivot].back()){
                                for(size_t ig=0; ig<grade_lists[pivot].size(); ig++){
                                    grade_t g = working_column.grade.join(grade_lists[pivot][ig]);
                                    if(visited_grades.find(g) == visited_grades.end()){
                                        grades.push(g);
                                        visited_grades.insert(g);
                                    }
                                }
                                grade_lists[pivot].push_back(working_column.grade);
                            }
                        }
                       
                    }else{
                        if(is_new){
                            gb_columns.pop_back();
                        }
                        working_column.syzygy.refresh();
                        Syz[i].push_back(signature_t(working_column.grade, -1));
                    }
                } else{
                    if(pivot != -1){
                        pivot_map[pivot] = column_index;
                        gb_columns[column_index].last_updated = iter_index;
                        if(is_new){
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, column_index));
                            if(grade_lists[pivot].size() == 0 || gb_columns[column_index].grade != grade_lists[pivot].back()){
                                for(size_t ig=0; ig<grade_lists[pivot].size(); ig++){
                                    grade_t g = gb_columns[column_index].grade.join(grade_lists[pivot][ig]);
                                    if(visited_grades.find(g) == visited_grades.end()){
                                        grades.push(g);
                                        visited_grades.insert(g);
                                    }
                                }
                                grade_lists[pivot].push_back(gb_columns[column_index].grade);
                            }
                        }
                    }
                }
            }
        }
    }
    std::cout << "Finished computing Groebner bases." << std::endl;
    return gb_columns;
}

Matrix buchberger(Matrix& columns){
    std::cout << "Computing GB for the image using Buchbergers algorithm." << std::endl;
    size_t max_pivot = 0;
    for(size_t i=0; i<columns.size(); i++){
        if(columns[i].get_pivot_index() > max_pivot){
            max_pivot = columns[i].get_pivot_index();
        }
    }
    std::vector<std::vector<size_t>> pivot_map(max_pivot+1, std::vector<size_t>());
    Matrix G = columns;
    std::queue<std::pair<size_t, size_t>> M;
    for(size_t i=0; i<G.size(); i++){
        for(size_t j=0; j<i; j++){
            if(G[i].get_pivot_index() == G[j].get_pivot_index()){
                M.push(std::pair<size_t, size_t>(i, j));
            }
        }
    }
    while(M.size()>0){
        std::pair<size_t, size_t> p = M.front();
        M.pop();
        SignatureColumn c = G[p.first];
        c.plus(G[p.second]);
        index_t pivot = c.get_pivot_index();
        while(pivot != -1){
            bool has_reduced = false;
            for(size_t i=0; i<pivot_map[pivot].size(); i++){
                if(G[pivot_map[pivot][i]].grade.leq_poset(c.grade)){
                    c.plus(G[pivot_map[pivot][i]]);
                    c.refresh();
                    c.syzygy.refresh();
                    has_reduced = true;
                    break;
                }
            }
            if(has_reduced){
                pivot = c.get_pivot_index();
            }else{
                
                for(size_t i=0; i<G.size(); i++){
                    if(G[i].get_pivot_index() == pivot){
                        M.push(std::pair<size_t, size_t>(i, G.size()));
                    }
                }
                G.push_back(c);
                pivot_map[pivot].push_back(G.size()-1);
                break;
            }
        }
        
    }
    std::cout << "Finished computing a GB of size: " << G.size() << std::endl;
    return G;
}

Matrix computeKernel_2p(std::vector<SignatureColumn>& columns){
    /*
     The main function computing a Groebner basis for the kernel of the bigraded map described by the list of columns 'columns'.
     
     Arguments:
     columns {std::vector<SignatureColumn>} -- columns describing the matrix of a map between two free multigraded momdules.
     
     Returns:
     std::vector<SignatureColumn> -- a list of vectors decribing a minimal Groebner basis for the image of the map.
     std::vector<SignatureColumn> -- a list of vectors describing a minimal Groebner basis for the kernel of the map.
     */
    
    std::cout << "Starting to compute Groebner bases..." << std::endl;
    
    /* Sort columns colexicographically */
    sort(columns.begin(), columns.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });
    // The sorted columns should agree with the columns sorted by index of signature
    hash_map<size_t, size_t> index_map_low, index_map_high;
    for(size_t i=0; i<columns.size(); i++){
        columns[i].signature_index = i;
        index_map_high[columns[i].grade[columns[i].grade.size()-1]] = i;
        if(index_map_low.find(columns[i].grade[columns[i].grade.size()-1]) == index_map_low.end()){
            index_map_low[columns[i].grade[columns[i].grade.size()-1]] = i;
        }
    }
    
    /* Compute index set iterator */
    std::vector<std::vector<index_t>> grade_base_set = get_grade_base_set<SignatureColumn>(columns);
    int n_total_grades = 1;
    for(auto& grade_list : grade_base_set){
        n_total_grades *= grade_list.size();
    }
    Iterator_lex grade_iterator = Iterator_lex(grade_base_set);
    
    std::vector<SyzColumn> syzygies;
    std::vector<std::vector<signature_t>> Syz;
    
    for(size_t i=0; i<columns.size(); i++){
        Syz.push_back(std::vector<signature_t>());
    }
    
    /* Main algorithm that iterates through the index set */
    index_t pivot;
    std::vector<index_t> pivot_map(columns.size(), -1);
    
    while(grade_iterator.has_next()){
        grade_t& v = grade_iterator.next();
        
        /* Initialize Macaulay matrix */
        for( size_t i=index_map_low[v[v.size()-1]]; i<=index_map_high[v[v.size()-1]] && columns[i].grade[0] <= v[0]; i++ ){
            if(columns[i].size()>0){
                pivot = columns[i].get_pivot_index();
                if(pivot_map[pivot] > -1 && pivot_map[pivot] < i){
                    SignatureColumn& working_column = columns[i];
                    while(pivot != -1 && pivot_map[pivot] > -1 && pivot_map[pivot] < i){
                        working_column.plus(columns[pivot_map[pivot]]);
                        pivot = working_column.get_pivot().get_index();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        pivot_map[pivot] = i;
                    }else{
                        working_column.syzygy.refresh();
                        syzygies.push_back(working_column.syzygy);
                        Syz[i].push_back(signature_t(working_column.grade, syzygies.size()-1));
                    }
                } else{
                    if(pivot != -1){
                        pivot_map[pivot] = i;
                    }
                }
            }
        }
    }
    
    Matrix syzygies_output;
    syzygies_output.reserve(syzygies.size());
    for(size_t i=0; i<Syz.size(); i++){
        for(size_t j=0; j<Syz[i].size(); j++){
            syzygies_output.push_back(SignatureColumn(Syz[i][j].get_grade(), syzygies_output.size(), syzygies[Syz[i][j].get_index()]));
        }
    }
    std::cout << "Finished computing Groebner bases." << std::endl;
    return syzygies_output;
}

/* Presentations */

hash_map<size_t, size_t> compute_local_pairs(Matrix& columns, hash_map<size_t, grade_t>& row_grades){
    /** Computes the local positive and negative pairs of columns and pivots.
     
     Arguments:
     columns {Matrix} -- the matrix with columns to be labeled local positive, negative or global.
     row_grades {std::vector<grade_t>} -- a list of the multigrade of each row.
     
     Returns:
     hash_map<size_t, size_t> -- a map that sends a local positive row to a local negative column.
     
     */
    std::cout << "Starting to compute local pairs...";
    hash_map<size_t, size_t> pairs;
    for(size_t column_index=0; column_index < columns.size(); column_index++){
        if(columns[column_index].local != 1 && columns[column_index].get_pivot().get_index() != -1 && row_grades[columns[column_index].get_pivot().get_index()] == columns[column_index].grade){
            
            if(pairs.find(columns[column_index].get_pivot().get_index()) != pairs.end()){
                SignatureColumn working_column(columns[column_index].grade, -1);
                size_t column_to_add = column_index;
                while(true){
                    working_column.plus(columns[column_to_add]);
                    index_t pivot = working_column.get_pivot().get_index();
                    if(pivot != -1 && row_grades[pivot] == working_column.grade){
                        if (pairs.find(pivot) != pairs.end()) {
                            column_to_add = pairs[pivot];
                        } else {
                            pairs[pivot] = column_index;
                            columns[column_index].local = -1;
                            break;
                        }
                    } else{
                        break;
                    }
                }
            }else{
                pairs[columns[column_index].get_pivot().get_index()] = column_index;
                columns[column_index].local = -1;
            }
        }
    }
    std::cout << "Finished computing local pairs.";
    return pairs;
}

Matrix compute_global_columns(Matrix& columns, hash_map<size_t, size_t>& positive_pairs, std::set<size_t>& negative_rows){
    /** Performs the presentation minimization step by removing local positive and negative columns.
     
     Arguments:
     columns {Matrix} -- a matrix with columns labeled local positive, negative or global.
     positive_pairs {hash_map<size_t, size_t>} -- a map sending local positive rows to local negative columns.
     negative_rows {std::set<size_t>} -- a set containing the indices of local negative rows.
     
     Returns:
     Matrix -- a minimized presentation matrix.
     
     */
    std::cout << "Starting to compute global columns...";
    
    Matrix global_columns;
    for(size_t index_column_to_reduce = 0; index_column_to_reduce<columns.size();index_column_to_reduce++) {
        if(columns[index_column_to_reduce].local != 0)
            continue;
        SignatureColumn working_boundary(columns[index_column_to_reduce].grade, index_column_to_reduce);
        SignatureColumn global_column(columns[index_column_to_reduce].grade, index_column_to_reduce);
        working_boundary.plus(columns[index_column_to_reduce]);
        while(true) {
            index_t pivot = working_boundary.get_pivot().get_index();
            if (pivot != -1) {
                if (negative_rows.find(pivot) != negative_rows.end()) {
                    working_boundary.pop_pivot();
                }else if(positive_pairs.find(pivot) != positive_pairs.end()) {
                    working_boundary.plus(columns[positive_pairs[pivot]]);
                }else{
                    global_column.push(working_boundary.pop_pivot());
                }
            }else{
                break;
            }
        }
        if(global_column.get_pivot().get_index() != -1){
            global_columns.push_back(global_column);
        }
    }
    std::cout << "Finished computing global columns of size: " << global_columns.size();
    return global_columns;
}


Matrix compute_syzygy_module(Matrix groebner_basis){
    /** Computes the syzygy module of a Groebner basis using Schreyer's algorithm.
     
    */
    
    for(size_t i=0; i<groebner_basis.size(); i++){
        groebner_basis[i].signature_index = i;
    }
    
    hash_map<size_t, std::vector<SignatureColumn>> pivot_partition;
    for(size_t i=0; i<groebner_basis.size(); i++){
        if(pivot_partition.find(groebner_basis[i].get_pivot().get_index()) == pivot_partition.end()){
            pivot_partition[groebner_basis[i].get_pivot().get_index()] = std::vector<SignatureColumn>();
        }
        pivot_partition[groebner_basis[i].get_pivot().get_index()].push_back(groebner_basis[i]);
    }
    Matrix syzygies;
    for(auto& entry : pivot_partition){
        std::vector<SignatureColumn>& columns = entry.second;
        //Matrix M;
        //for(auto& column : columns){
        //    M.push_back(column);
        //}
        //std::cout << "Partition for pivot " << entry.first << "\n";
        //M.print();
        for(size_t i=1; i<columns.size(); i++){
            std::vector<grade_t> minimal_elements_tmp;
            std::vector<size_t> minimal_indices_tmp;
            for(size_t j=0; j<i; j++){
                grade_t m_ji = columns[j].grade.m_ji(columns[i].grade);
                //spdlog::info("m_ji...");
                //columns[j].grade.print();
                //columns[i].grade.print();
                //m_ji.print();
                bool is_minimal = true;
                for(auto& el : minimal_elements_tmp){
                    if(el.leq_poset(m_ji)){
                        is_minimal = false;
                        break;
                    }
                }
                //std::cout << "Is minimal: " << is_minimal << "\n";
                if(is_minimal){
                    minimal_elements_tmp.push_back(m_ji);
                    minimal_indices_tmp.push_back(j);
                }
            }
            std::vector<grade_t> minimal_elements;
            std::vector<size_t> minimal_indices;
            for(size_t j=0; j<minimal_elements_tmp.size(); j++){
                bool is_minimal = true;
                for(size_t k=0; k<minimal_elements_tmp.size(); k++){
                    if(k!=j && minimal_elements_tmp[k].leq_poset(minimal_elements_tmp[j])){
                        is_minimal = false;
                        break;
                    }
                }
                if(is_minimal){
                    minimal_elements.push_back(minimal_elements_tmp[j]);
                    minimal_indices.push_back(minimal_indices_tmp[j]);
                }
            }
            /*spdlog::info("Minimal elements...");
            for(auto& min_el : minimal_elements){
                min_el.print();
            }*/
            for(size_t j=0; j<minimal_elements.size(); j++){
                SignatureColumn column_to_reduce(columns[i]);
                column_to_reduce.plus(columns[minimal_indices[j]]);
                SignatureColumn syzygy(column_to_reduce.grade, -1);
                index_t pivot = column_to_reduce.get_pivot().get_index();
                while(pivot != -1){
                    if(pivot_partition.find(pivot) == pivot_partition.end()){
                        throw "Unkown pivot.";
                    }
                    bool has_reduced = false;
                    for(auto& column : pivot_partition[pivot]){
                        if(column.grade.leq_poset(column_to_reduce.grade)){
                            column_to_reduce.plus(column);
                            syzygy.push(column_entry_t(1, column.signature_index));
                            has_reduced = true;
                            break;
                        }
                    }
                    if(!has_reduced){
                        throw "Could not reduce column.";
                    }
                    pivot = column_to_reduce.get_pivot().get_index();
                }
                syzygy.push(column_entry_t(1, columns[i].signature_index));
                syzygy.push(column_entry_t(1, columns[minimal_indices[j]].signature_index));
                if(syzygy.get_pivot().get_index() != -1){
                    syzygies.push_back(syzygy);
                }
                
            }
        }
    }
    return syzygies;
}

Matrix compute_minimal_generating_set(Matrix& generators){
    /** Computes a minimal generating set for the module described by the columns in 'generators'.
     
     Arguments:
     generators {Matrix} -- a matrix whose columns a generators of a module.
     
     Returns:
     Matrix -- a subset of the columns of 'generators' that constitute a minimal set of generators for the module.
     */
    
    /* Vectors to store the columns of the GBs */
    std::cout << "Starting to compute minimal generating set for columns of size: " << generators.size() << std::endl;
    
    /* Sort columns colexicographically */
    sort(generators.begin(), generators.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });
    
    Matrix gb_columns;
    
    /* Index maps to keep track of signatures */
    std::vector<std::vector<signature_t>> GB; //The columns of 'generators' with a non-empty signature list is a generator.
    std::vector<std::vector<signature_t>> Syz;
    
    for(size_t i=0; i<generators.size(); i++){
        GB.push_back(std::vector<signature_t>());
        Syz.push_back(std::vector<signature_t>());
        generators[i].signature_index = i;
    }
    
    hash_map<size_t, size_t> index_map_high;
    for(size_t i=0; i<generators.size(); i++){
        generators[i].signature_index = i;
        index_map_high[generators[i].grade[generators[i].grade.size()-1]] = i;
    }
    
    std::vector<grade_t> index_list;
    for(size_t column_index=0 ; column_index<generators.size(); column_index++){
        index_list.push_back(generators[column_index].grade);
    }
    sort(index_list.begin(), index_list.end()); // Sort grades lexicographically
    
    
    index_t max_pivot=0;
    for(auto& generator : generators){
        if(max_pivot < generator.get_pivot_index()){
            max_pivot = generator.get_pivot_index();
        }
    }
    
    std::vector<size_t> grade_hashes;
    grade_hashes.reserve(generators.size());
    GradeHasher grade_hasher;
    for(auto& c : generators){
        grade_hashes.push_back(grade_hasher(c.grade));
    }
    
    /* Main algorithm that iterates through the index set */
    int column_index;
    index_t pivot;
    std::vector<index_t> pivot_map(max_pivot+1, -1);
    int iter_index = 0;
    
    for(size_t index=0 ; index<index_list.size(); index++){
        if(index > 0 && index_list[index]==index_list[index-1]){
            continue;
        }
        grade_t v = index_list[index];
        iter_index++;
        size_t grade_hash = grade_hasher(v);
        
        
        /* Initialize Macaulay matrix */
        for( size_t i=0; i<=index_map_high[v[v.size()-1]]; i++ ){
            column_index = -1;
            bool is_new=false;
            if(GB[i].size() > 0){
                bool in_syz = false;
                if(Syz[i].size()>0){
                    for(size_t j=Syz[i].size()-1; j < Syz[i].size(); j--){
                        if((Syz[i][j].get_grade()).leq_poset(v)){
                            in_syz = true;
                            break;
                        }
                    }
                }
                if(!in_syz){
                    for(size_t j=GB[i].size()-1; j < GB[i].size(); j--){
                        if(GB[i][j].get_grade().leq_poset(v)){
                            column_index = (int)GB[i][j].get_index();
                            break;
                        }
                    }
                }
            } else{
                if(grade_hash == grade_hashes[i] && generators[i].grade == v){
                    gb_columns.push_back(generators[i]);
                    column_index = (int)gb_columns.size()-1;
                    is_new=true;
                }
            }
            if(column_index > -1){
                pivot = gb_columns[column_index].get_pivot().get_index();
                if(pivot > -1 && pivot_map[pivot] > -1 && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                    SignatureColumn working_column = gb_columns[column_index];
                    while(pivot != -1 && pivot_map[pivot] > -1 && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                        working_column.plus(gb_columns[pivot_map[pivot]]);
                        pivot = working_column.get_pivot().get_index();
                    }
                    if(is_new){
                        gb_columns.pop_back();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        gb_columns.push_back(working_column);
                        GB[i].push_back(signature_t(working_column.grade, gb_columns.size()-1));
                        pivot_map[pivot] = gb_columns.size()-1;
                        gb_columns[gb_columns.size()-1].last_updated = iter_index;
                    }else{
                        working_column.syzygy.refresh();
                        Syz[i].push_back(signature_t(working_column.grade, -1));
                    }
                } else{
                    if(pivot != -1){
                        pivot_map[pivot] = column_index;
                        gb_columns[column_index].last_updated = iter_index;
                        if(is_new){
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, column_index));
                        }
                    }
                }
            }
        }
    }
    Matrix minimal_generating_set;
    for(size_t column_index=0 ; column_index<generators.size(); column_index++){
        if(GB[column_index].size()>0){
            minimal_generating_set.push_back(generators[column_index]);
        }
    }
    //get_mem_usage(virt_memory, res_memory);
    return minimal_generating_set;
}

std::pair<Matrix, Matrix> compute_minimal_generating_set2(Matrix generators){
    /** Computes a minimal generating set for the module described by the columns in 'generators'.
     
     Arguments:
     generators {Matrix} -- a matrix whose columns a generators of a module.
     
     Returns:
     Matrix -- a subset of the columns of 'generators' that constitute a minimal set of generators for the module.
     */
    
    /* Vectors to store the columns of the GBs */
    std::cout << "Starting to compute minimal generating set for columns of size: " << generators.size()  << std::endl;
    
    for(size_t i=0; i<generators.size(); i++){
        generators[i].signature_index = i;
    }
    
    /* Sort columns colexicographically */
    sort(generators.begin(), generators.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });
    
    hash_map<size_t, size_t> reorder_map;
    for(size_t i=0; i<generators.size(); i++){
        reorder_map[i] = generators[i].signature_index;
    }
    
    Matrix gb_columns;
    std::vector<SyzColumn> syzygies;
    
    /* Index maps to keep track of signatures */
    std::vector<std::vector<signature_t>> GB; //The columns of 'generators' with a non-empty signature list is a generator.
    std::vector<std::vector<signature_t>> Syz;
    
    for(size_t i=0; i<generators.size(); i++){
        GB.push_back(std::vector<signature_t>());
        Syz.push_back(std::vector<signature_t>());
        generators[i].signature_index = i;
        generators[i].syzygy = SyzColumn();
        generators[i].syzygy.push(column_entry_t(1, i));
    }
    
    hash_map<size_t, size_t> index_map_high;
    for(size_t i=0; i<generators.size(); i++){
        index_map_high[generators[i].grade[generators[i].grade.size()-1]] = i;
    }
    
    std::vector<grade_t> index_list;
    for(size_t column_index=0 ; column_index<generators.size(); column_index++){
        index_list.push_back(generators[column_index].grade);
    }
    sort(index_list.begin(), index_list.end()); // Sort grades lexicographically
    
    index_t max_pivot=0;
    for(auto& generator : generators){
        if(max_pivot < generator.get_pivot().get_index()){
            max_pivot = generator.get_pivot().get_index();
        }
    }
    
    std::vector<size_t> grade_hashes;
    grade_hashes.reserve(generators.size());
    GradeHasher grade_hasher;
    for(auto& c : generators){
        grade_hashes.push_back(grade_hasher(c.grade));
    }
    
    /* Main algorithm that iterates through the index set */
    int column_index;
    index_t pivot;
    
    int iter_index = 0;
    
    for(size_t index=0 ; index<index_list.size(); index++){
       // if(index > 0 && index_list[index]==index_list[index-1]){
       //     continue;
       // }
        grade_t& v = index_list[index];
        iter_index++;
        std::vector<index_t> pivot_map(max_pivot+1, -1);
        size_t grade_hash = grade_hasher(v);
        /* Initialize Macaulay matrix */
        for( size_t i=0; i<generators.size(); i++ ){
            column_index = -1;
            bool is_new=false;
            if(GB[i].size() > 0){
                bool in_syz = false;
                if(Syz[i].size()>0){
                    for(size_t j=Syz[i].size()-1; j < Syz[i].size(); j--){
                        if((Syz[i][j].get_grade()).leq_poset(v)){
                            in_syz = true;
                            break;
                        }
                    }
                }
                if(!in_syz){
                    for(size_t j=GB[i].size()-1; j < GB[i].size(); j--){
                        if(GB[i][j].get_grade().leq_poset(v)){
                            column_index = (int)GB[i][j].get_index();
                            break;
                        }
                    }
                }
            } else{
                if(grade_hash == grade_hashes[i] && generators[i].grade == v){
                    gb_columns.push_back(generators[i]);
                    column_index = (int)gb_columns.size()-1;
                    is_new=true;
                }
            }
            if(column_index > -1){
                pivot = gb_columns[column_index].get_pivot().get_index();
                if(pivot > -1 && pivot_map[pivot] > -1){
                    SignatureColumn working_column(gb_columns[column_index]);
                    SyzColumn syz;
                    while(pivot != -1 && pivot_map[pivot] > -1){
                        working_column.plus(gb_columns[pivot_map[pivot]]);
                        syz.plus(gb_columns[pivot_map[pivot]].syzygy);
                        pivot = working_column.get_pivot().get_index();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        gb_columns.push_back(working_column);
                        GB[i].push_back(signature_t(working_column.grade, gb_columns.size()-1));
                        pivot_map[pivot] = gb_columns.size()-1;
                        gb_columns[gb_columns.size()-1].last_updated = iter_index;
                    }else{
                        working_column.syzygy.refresh();
                        syz.refresh();
                        syzygies.push_back(syz);
                        Syz[i].push_back(signature_t(working_column.grade, syzygies.size()-1));
                    }
                } else{
                    if(pivot != -1){
                        pivot_map[pivot] = column_index;
                        gb_columns[column_index].last_updated = iter_index;
                        if(is_new){
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, column_index));
                        }
                    }
                }
            }
        }
    }
    Matrix minimal_generating_set;
    Matrix change_of_basis_map;
    hash_map<size_t, size_t> reindex_map;
    for(size_t column_index=0 ; column_index<generators.size(); column_index++){
        if(GB[column_index].size()>0){
            reindex_map[column_index] = minimal_generating_set.size();
            minimal_generating_set.push_back(generators[column_index]);
            SignatureColumn column(generators[column_index].grade, reorder_map[column_index]);
            column.push(column_entry_t(1, reindex_map[column_index]));
            change_of_basis_map.push_back(column);
        }else{
            SignatureColumn column(generators[column_index].grade, reorder_map[column_index]);
            index_t pivot = syzygies[Syz[column_index][0].get_index()].get_pivot_index();
            while(pivot != -1){
                column.push(column_entry_t(1, reindex_map[pivot]));
                syzygies[Syz[column_index][0].get_index()].pop_pivot();
                pivot = syzygies[Syz[column_index][0].get_index()].get_pivot_index();
            }
            change_of_basis_map.push_back(column);
        }
    }
    
    sort(change_of_basis_map.begin(), change_of_basis_map.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.signature_index < rhs.signature_index;
         });
    get_mem_usage(virt_memory, res_memory);
    std::cout << "Finished computing minimal generating set of size: " << minimal_generating_set.size();
    return std::pair<Matrix, Matrix>(minimal_generating_set, change_of_basis_map);
}


bool cmpID(const signature_t& lhs, const signature_t& rhs)
{
    return lhs.first.lt_colex(rhs.first);
}

void insert_sorted( std::vector<signature_t>& cont, signature_t value ) {
    std::vector<signature_t>::iterator it = std::upper_bound( cont.begin(), cont.end(), value, cmpID); // find proper position in descending order
    cont.insert( it, value ); // insert before iterator it
}

std::pair<Matrix, hash_map<size_t, grade_t>> computePresentationDeg_imopt(Matrix& image_columns, Matrix& columns, bool debug=true){
    /*
     The main function computing a Groebner basis for the image and kernel of the map described by the list of columns 'columns'.
     
     Arguments:
     image_columns {std::vector<SignatureColumn>} -- a minimal generating set of the image map.
     columns {std::vector<SignatureColumn>} -- columns describing the matrix of a map between two free multigraded modules.
     
     Returns:
     std::vector<SignatureColumn> -- a list of vectors decribing a minimal Groebner basis for the image of the map.
     std::vector<SignatureColumn> -- a list of vectors describing a minimal Groebner basis for the kernel of the map.
     */
    
    std::cout << "Starting to presentation degree by degree..."  << std::endl;
    
    /* Sort columns colexicographically */
    sort(columns.begin(), columns.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });
    sort(image_columns.begin(), image_columns.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade<rhs.grade;
         });
    // The sorted columns should agree with the columns sorted by index of signature
    hash_map<size_t, size_t> index_map_high;
    for(size_t i=0; i<columns.size(); i++){
        columns[i].signature_index = i;
        columns[i].syzygy = SyzColumn();
        columns[i].syzygy.push(column_entry_t(1, i));
        index_map_high[columns[i].grade[columns[i].grade.size()-1]] = i;
    }
    
    /* Compute index set iterator */
    std::priority_queue<grade_t, std::vector<grade_t>, std::greater<grade_t>> grades;
    std::vector<std::vector<grade_t>> grade_lists;
    std::vector<std::vector<grade_t>> grade_listsH;
    
    /* Vectors to store the columns of the GBs */
    Matrix gb_columns, gb_columnsH, syzygiesH;
    
    
    /* Index maps to keep track of signatures */
    std::vector<std::vector<signature_t>> GB, GBH;
    std::vector<std::vector<signature_t>> Syz, SyzH;
    
    for(size_t i=0; i<columns.size(); i++){
        GB.push_back(std::vector<signature_t>());
        Syz.push_back(std::vector<signature_t>());
    }
    
    /* Main algorithm that iterates through the index set */
    int column_index;
    index_t max_pivot=0;
    for(auto& column : columns){
        if(max_pivot < column.get_pivot_index()){
            max_pivot = column.get_pivot_index();
        }
    }
    for(size_t i=0; i<=max_pivot; i++){
        grade_lists.push_back(std::vector<grade_t>());
    }
    for(auto& column : columns){
        grades.push(column.grade);
        grade_listsH.push_back(std::vector<grade_t>());
    }
    for(auto& column : image_columns){
        grades.push(column.grade);
    }
    index_t pivot;
    
    std::vector<index_t> pivot_map(max_pivot+1, -1);
    std::vector<index_t> pivot_mapH(columns.size(), -1);
    
    std::vector<size_t> grade_hashes;
    grade_hashes.reserve(columns.size());
    GradeHasher grade_hasher;
    for(auto& c : columns){
        grade_hashes.push_back(grade_hasher(c.grade));
    }
    
    std::vector<signature_t> reorder_Z_map;
    
    std::set<size_t> syz_pivots;
    
    size_t image_index=0;
    
    
    int iter_index = 0;
    while(!grades.empty()){
        grade_t v = grades.top();
        while(v == grades.top()){
            grades.pop();
        }
    
        iter_index++;
        size_t grade_hash = grade_hasher(v);
        
        // B^h_z
        for( auto& signature : reorder_Z_map ){
            size_t i = signature.get_index();
            column_index = -1;
            //bool in_syz = false;
            //size_t syz_size = SyzH[i].size();
            if(SyzH[i].size()==0 && GBH[i][0].get_grade().leq_poset(v)){
                size_t gb_size = GBH[i].size();
                column_index = (int)GBH[i][0].get_index();
                for(size_t j=gb_size-1; j >= 1; j--){
                    if(GBH[i][j].get_grade().leq_poset(v)){
                        column_index = (int)GBH[i][j].get_index();
                        break;
                    }
                }
                /*for(size_t j=syz_size-1; j < syz_size; j--){
                    if((SyzH[i][j].get_grade()).leq_poset(v)){
                        in_syz = true;
                        break;
                    }
                }*/
            }
            
            /*if(!in_syz){
                size_t gb_size = GBH[i].size();
                for(size_t j=gb_size-1; j < gb_size; j--){
                    if(GBH[i][j].get_grade().leq_poset(v)){
                        column_index = (int)GBH[i][j].get_index();
                        break;
                    }
                }
            }*/
            
            
            if(column_index > -1){
                pivot = gb_columnsH[column_index].get_pivot_index();
                if(pivot > -1 && pivot_mapH[pivot] > -1 && gb_columnsH[pivot_mapH[pivot]].last_updated == iter_index){
                    SignatureColumn working_column = gb_columnsH[column_index];
                    while(pivot != -1 && pivot_mapH[pivot] > -1 && gb_columnsH[pivot_mapH[pivot]].last_updated == iter_index){
                        working_column.plus(gb_columnsH[pivot_mapH[pivot]]);
                        pivot = working_column.get_pivot_index();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        grade_t grade = working_column.grade;
                        /*if(GBH[i].size()==1){
                            gb_columnsH.push_back(working_column);
                            GBH[i].push_back(signature_t(grade, gb_columnsH.size()-1));
                        }else{
                            if(pivot_mapH[gb_columnsH[GBH[i][1].get_index()].get_pivot_index()] == GBH[i][1].get_index()){
                                pivot_mapH[gb_columnsH[GBH[i][1].get_index()].get_pivot_index()] = -1;
                            }
                            gb_columnsH[GBH[i][1].get_index()] = working_column;
                            GBH[i][1].first = grade;
                        }*/
                        gb_columnsH.push_back(working_column);
                        GBH[i].push_back(signature_t(grade, gb_columnsH.size()-1));
                        pivot_mapH[pivot] = gb_columnsH.size()-1;
                        gb_columnsH[gb_columnsH.size()-1].last_updated = iter_index;
                        if(grade == v){
                            if(grade_listsH[pivot].size() == 0 || grade != grade_listsH[pivot].back()){
                                std::vector<grade_t> minimal_elements_tmp;
                                for(size_t j=0; j<grade_listsH[pivot].size(); j++){
                                    grade_t m_ji = grade_listsH[pivot][j].m_ji(grade);
                                    bool is_minimal = true;
                                    for(auto& el : minimal_elements_tmp){
                                        if(el.leq_poset_m(m_ji)){
                                            is_minimal = false;
                                            break;
                                        }
                                    }
                                    if(is_minimal){
                                        minimal_elements_tmp.push_back(m_ji);
                                        grades.push(grade.join(grade_listsH[pivot][j]));
                                    }
                                }
                                grade_listsH[pivot].push_back(grade);
                            }
                        }
                    }else{
                        working_column.syzygy.refresh();
                       // if(syz_pivots.find(i) == syz_pivots.end()){
                            syzygiesH.push_back(SignatureColumn(working_column.grade, syzygiesH.size(), working_column.syzygy));
                            SyzH[i].push_back(signature_t(working_column.grade, syzygiesH.size()-1));
                            syzygiesH[syzygiesH.size()-1].last_updated = iter_index;
                      /*  }else{
                            SyzH[i].push_back(signature_t(working_column.grade, syzygiesH.size()-1));
                        }*/
                    }
                } else{
                    if(pivot != -1){
                        pivot_mapH[pivot] = column_index;
                        gb_columnsH[column_index].last_updated = iter_index;
                    }
                }
            }
        }
        
        /* Reduce image columns */
        while(image_index < image_columns.size() && image_columns[image_index].grade == v){
            SignatureColumn working_column(image_columns[image_index]);
            working_column.syzygy = SyzColumn();
            pivot = working_column.get_pivot_index();
            while(pivot != -1 && pivot_mapH[pivot] != -1 && gb_columnsH[pivot_mapH[pivot]].last_updated == iter_index){
                working_column.plus(gb_columnsH[pivot_mapH[pivot]]);
                pivot = working_column.get_pivot_index();
            }
            working_column.syzygy.refresh();
            if(pivot != -1){
                //throw "Found non-reduced image column";
                working_column.refresh();
                working_column.signature_index = GBH.size();
                GBH.push_back(std::vector<signature_t>());
                SyzH.push_back(std::vector<signature_t>());
                gb_columnsH.push_back(working_column);
                GBH[GBH.size()-1].push_back(signature_t(working_column.grade, gb_columnsH.size()-1));
                insert_sorted(reorder_Z_map, signature_t(working_column.grade, GBH.size()-1));
                pivot_mapH[pivot] = gb_columnsH.size()-1;
                gb_columnsH[gb_columnsH.size()-1].last_updated = iter_index;
                if(grade_listsH[pivot].size() == 0 || working_column.grade != grade_listsH[pivot].back()){
                    std::vector<grade_t> minimal_elements_tmp;
                    for(size_t j=0; j<grade_listsH[pivot].size(); j++){
                        grade_t m_ji = grade_listsH[pivot][j].m_ji(working_column.grade);
                        bool is_minimal = true;
                        for(auto& el : minimal_elements_tmp){
                            if(el.leq_poset_m(m_ji)){
                                is_minimal = false;
                                break;
                            }
                        }
                        if(is_minimal){
                            minimal_elements_tmp.push_back(m_ji);
                            grades.push(working_column.grade.join(grade_listsH[pivot][j]));
                        }
                    }
                    grade_listsH[pivot].push_back(working_column.grade);
                }
                //working_column.syzygy = SyzColumn();
                //working_column.syzygy.push(column_entry_t(1, gb_columnsH.size()-1));
                //syzygiesH.push_back(SignatureColumn(working_column.grade, syzygiesH.size(), working_column.syzygy));
            }else{
                syzygiesH.push_back(SignatureColumn(working_column.grade, syzygiesH.size(), working_column.syzygy));
                //syz_pivots.insert(working_column.syzygy.get_pivot_index());
            }
            image_index++;
        }
        
        /* Initialize Macaulay matrix */
        size_t& index_bound = index_map_high[v[v.size()-1]];
        for( size_t i=0; i<=index_bound; i++ ){
            column_index = -1;
            bool is_new=false;
            if(GB[i].size() > 0){
                if(pivot_mapH[i] == -1 || gb_columnsH[pivot_mapH[i]].last_updated != iter_index){
                    for(size_t j=GB[i].size()-1; j < GB[i].size(); j--){
                        if(GB[i][j].get_grade().leq_poset(v)){
                            column_index = (int)GB[i][j].get_index();
                            break;
                        }
                    }
                }
            } else{
                if(pivot_mapH[i] == -1 && grade_hash == grade_hashes[i] && columns[i].grade == v){
                    gb_columns.push_back(columns[i]);
                    column_index = (int)gb_columns.size()-1;
                    is_new=true;
                }
            }
            if(column_index > -1){
                pivot = gb_columns[column_index].get_pivot_index();
                if(pivot != -1 && pivot_map[pivot] != -1 && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                    SignatureColumn working_column = gb_columns[column_index];
                    while(pivot != -1 && pivot_map[pivot] != -1 && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                        working_column.plus(gb_columns[pivot_map[pivot]]);
                        pivot = working_column.get_pivot_index();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        grade_t grade = working_column.grade;
                        if(is_new){
                            gb_columns[column_index] = working_column;
                            GB[i].push_back(signature_t(grade, column_index));
                            pivot_map[pivot] = column_index;
                            gb_columns[column_index].last_updated = iter_index;
                        }else{
                            /*if(GB[i].size()==1){
                                gb_columns.push_back(working_column);
                                GB[i].push_back(signature_t(grade, gb_columns.size()-1));
                            }else{
                                if(pivot_map[gb_columns[GB[i][1].get_index()].get_pivot_index()] == GB[i][1].get_index()){
                                    pivot_map[gb_columns[GB[i][1].get_index()].get_pivot_index()] = -1;
                                }
                                gb_columns[GB[i][1].get_index()] = working_column;
                                GB[i][1].first = grade;
                            }*/
                            gb_columns.push_back(working_column);
                            GB[i].push_back(signature_t(grade, gb_columns.size()-1));
                            pivot_map[pivot] = gb_columns.size()-1;
                            gb_columns[gb_columns.size()-1].last_updated = iter_index;
                        }
                        if(grade == v){
                            if(grade_lists[pivot].size() == 0 || grade != grade_lists[pivot].back()){
                                std::vector<grade_t> minimal_elements_tmp;
                                for(size_t j=0; j<grade_lists[pivot].size(); j++){
                                    grade_t m_ji = grade_lists[pivot][j].m_ji(grade);
                                    bool is_minimal = true;
                                    for(auto& el : minimal_elements_tmp){
                                        if(el.leq_poset_m(m_ji)){
                                            is_minimal = false;
                                            break;
                                        }
                                    }
                                    if(is_minimal){
                                        minimal_elements_tmp.push_back(m_ji);
                                        grades.push(grade.join(grade_lists[pivot][j]));
                                    }
                                }
                                grade_lists[pivot].push_back(grade);
                            }
                        }
                    }else{
                        if(is_new){
                            gb_columns.pop_back();
                        }
                        working_column.syzygy.refresh();
                        SignatureColumn syz(working_column.grade, GBH.size(), working_column.syzygy);
                        syz.syzygy.push(column_entry_t(1, GBH.size()));
                        GBH.push_back(std::vector<signature_t>());
                        SyzH.push_back(std::vector<signature_t>());
                        gb_columnsH.push_back(syz);
                        GBH[GBH.size()-1].push_back(signature_t(syz.grade, gb_columnsH.size()-1));
                        insert_sorted(reorder_Z_map, signature_t(working_column.grade, GBH.size()-1));
                        pivot_mapH[i] = gb_columnsH.size()-1;
                        gb_columnsH[gb_columnsH.size()-1].last_updated = iter_index;
                        if(working_column.grade == v){
                            if(grade_listsH[i].size() == 0 || working_column.grade != grade_listsH[i].back()){
                                std::vector<grade_t> minimal_elements_tmp;
                                for(size_t j=0; j<grade_listsH[i].size(); j++){
                                    grade_t m_ji = grade_listsH[i][j].m_ji(working_column.grade);
                                    bool is_minimal = true;
                                    for(auto& el : minimal_elements_tmp){
                                        if(el.leq_poset_m(m_ji)){
                                            is_minimal = false;
                                            break;
                                        }
                                    }
                                    if(is_minimal){
                                        minimal_elements_tmp.push_back(m_ji);
                                        grades.push(working_column.grade.join(grade_listsH[i][j]));
                                    }
                                }
                                grade_listsH[i].push_back(working_column.grade);
                            }
                        }
                    }
                } else{
                    if(pivot != -1){
                        if(is_new){
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, column_index));
                            pivot_map[pivot] = column_index;
                            gb_columns[column_index].last_updated = iter_index;
                            if(grade_lists[pivot].size() == 0 || gb_columns[column_index].grade != grade_lists[pivot].back()){
                                std::vector<grade_t> minimal_elements_tmp;
                                for(size_t j=0; j<grade_lists[pivot].size(); j++){
                                    grade_t m_ji = grade_lists[pivot][j].m_ji(gb_columns[column_index].grade);
                                    bool is_minimal = true;
                                    for(auto& el : minimal_elements_tmp){
                                        if(el.leq_poset_m(m_ji)){
                                            is_minimal = false;
                                            break;
                                        }
                                    }
                                    if(is_minimal){
                                        minimal_elements_tmp.push_back(m_ji);
                                        grades.push(gb_columns[column_index].grade.join(grade_lists[pivot][j]));
                                    }
                                }
                                grade_lists[pivot].push_back(gb_columns[column_index].grade);
                            }
                        }else{
                            pivot_map[pivot] = column_index;
                            gb_columns[column_index].last_updated = iter_index;
                        }
                    } else {
                        if(is_new){
                            SignatureColumn working_column = gb_columns[column_index];
                            SyzColumn syz_column = SyzColumn();
                            syz_column.push(column_entry_t(1, column_index));
                            SignatureColumn syz(working_column.grade, GBH.size(), syz_column);
                            syz.syzygy.push(column_entry_t(1, GBH.size()));
                            GBH.push_back(std::vector<signature_t>());
                            SyzH.push_back(std::vector<signature_t>());
                            gb_columnsH.push_back(syz);
                            GBH[GBH.size()-1].push_back(signature_t(syz.grade, gb_columnsH.size()-1));
                            insert_sorted(reorder_Z_map, signature_t(working_column.grade, GBH.size()-1));
                            pivot_mapH[i] = gb_columnsH.size()-1;
                            gb_columnsH[gb_columnsH.size()-1].last_updated = iter_index;
                            if(working_column.grade == v){
                                if(grade_listsH[i].size() == 0 || working_column.grade != grade_listsH[i].back()){
                                    std::vector<grade_t> minimal_elements_tmp;
                                    for(size_t j=0; j<grade_listsH[i].size(); j++){
                                        grade_t m_ji = grade_listsH[i][j].m_ji(working_column.grade);
                                        bool is_minimal = true;
                                        for(auto& el : minimal_elements_tmp){
                                            if(el.leq_poset_m(m_ji)){
                                                is_minimal = false;
                                                break;
                                            }
                                        }
                                        if(is_minimal){
                                            minimal_elements_tmp.push_back(m_ji);
                                            grades.push(working_column.grade.join(grade_listsH[i][j]));
                                        }
                                    }
                                    grade_listsH[i].push_back(working_column.grade);
                                }
                            }
                        
                        }
                    }
                }
            }
        }
    }
    
    hash_map<size_t, grade_t> row_grade_map;
    for(size_t j=0; j<gb_columnsH.size(); j++){
        row_grade_map[j] = grade_t(gb_columnsH[j].grade);
    }

    std::cout << "Nonminimal presentation of size "<< syzygiesH.size() << std::endl;
    
    
    std::vector<index_t> rows;
    for(auto& entry : row_grade_map){
        rows.push_back(entry.first);
    }
    sort(rows.begin(), rows.end());
    hash_map<size_t, size_t> row_index_map;
    for(size_t i=0; i<rows.size(); i++){
        row_index_map[rows[i]] = i;
    }

    syzygiesH = compute_minimal_generating_set(syzygiesH);

//get_mem_usage(virt_memory, res_memory);
    return std::pair<Matrix, hash_map<size_t, grade_t>>(syzygiesH, row_grade_map);
}

std::pair<Matrix, hash_map<size_t, grade_t>> compute_presentation_2p(Matrix& image_generators, Matrix& generating_set_kernel){
    /** Computes a minimal presentation of the ith homology module using generators of the kernel and image of
     the boundary matrices \Delta_{i-1} and \Delta_i respectively.
     
     Arguments:
     image_generators {Matrix} -- a generating set of the image of the boundary map \Delta_i.
     kernel_generators {Matrix} -- a Groebner basis of the kernel of the boundary matrix \Delta_{i-1}.
     
     Returns:
     Matrix -- a matrix describing a minimal presentation of the ith homology module.
     
     */
    
    std::cout << "Starting to compute presentation...";
    Matrix generating_set_image = compute_minimal_generating_set(image_generators);
    
    //Sort image and kernel gens in colex and lex respectively
   sort(generating_set_kernel.begin(), generating_set_kernel.end(),[ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade) ;
         });
     sort(generating_set_image.begin(), generating_set_image.end(),[ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade < rhs.grade ;
         });
    
    generating_set_kernel.print();
    hash_map<size_t, grade_t> row_grade_map;
    for(size_t j=0; j<generating_set_kernel.size(); j++){
        row_grade_map[j] = grade_t(generating_set_kernel[j].grade);
        row_grade_map[j].print();
    }
    for(size_t i=0; i<generating_set_kernel.size(); i++){
        generating_set_kernel[i].syzygy = SyzColumn();
        generating_set_kernel[i].syzygy.push(column_entry_t(1, i));
    }
    
    std::cout << "Expressing the image columns in terms of the kernel...";
    
    hash_map<size_t, size_t> index_map_low, index_map_high;
    for(size_t i=0; i<generating_set_kernel.size(); i++){
        generating_set_kernel[i].signature_index = i;
        index_map_high[generating_set_kernel[i].grade.back()] = i;
        if(index_map_low.find(generating_set_kernel[i].grade.back()) == index_map_low.end()){
            index_map_low[generating_set_kernel[i].grade.back()] = i;
        }
    }
    
    
    /* Main algorithm that iterates through the index set */
    index_t max_pivot=0;
    for(auto& generator : generating_set_kernel){
        if(max_pivot < generator.get_pivot().get_index()){
            max_pivot = generator.get_pivot().get_index();
        }
    }
    index_t pivot;
    std::vector<index_t> pivot_map(max_pivot+1, -1);
    Matrix presentation_matrix;
    presentation_matrix.reserve(generating_set_image.size());
    for(auto& column : generating_set_image){
        grade_t& v = column.grade;
        
        /* Initialize Macaulay matrix */
        for( size_t i=0; i<=index_map_high[v.back()]; i++ ){
            if(generating_set_kernel[i].size()>0 && generating_set_kernel[i].grade[0] <= v[0]){
                pivot = generating_set_kernel[i].get_pivot().get_index();
                if(pivot_map[pivot] > -1 && pivot_map[pivot] < i){
                    SignatureColumn& working_column = generating_set_kernel[i];
                    while(pivot != -1 && pivot_map[pivot] > -1 && pivot_map[pivot] < i){
                        working_column.plus(generating_set_kernel[pivot_map[pivot]]);
                        pivot = working_column.get_pivot().get_index();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        pivot_map[pivot] = i;
                    }
                } else{
                    if(pivot != -1){
                        pivot_map[pivot] = i;
                    }
                }
            }
        }
        
        
        // Reduce the image column with the pivot set.
        column.syzygy = SyzColumn();
        while(true){
            pivot = column.get_pivot().get_index();
            if(pivot != -1){
                if(pivot_map[pivot] == -1){
                    std::cerr << "Cannot express image column in terms of kernel. Throwing exception...";
                    throw "Failed to express image column in terms of kernel generating set.";
                }
                column.plus(generating_set_kernel[pivot_map[pivot]]);
            }else{
                break;
            }
        }
        presentation_matrix.push_back(SignatureColumn(column.grade, column.signature_index, column.syzygy));
        presentation_matrix.back().refresh();
    }
    
    
    presentation_matrix.print();
    for(auto& column : presentation_matrix){
        column.grade.print();
    }
    
    hash_map<size_t, size_t> pairs = compute_local_pairs(presentation_matrix, row_grade_map);
    
    std::set<size_t> negative_columns;
    presentation_matrix = compute_global_columns(presentation_matrix, pairs, negative_columns);
    
    for(auto& entry : pairs){
        if(row_grade_map.find(entry.first) != row_grade_map.end()){
            row_grade_map.erase(entry.first);
        }
    }
    
    
    
    // Reindex the rows of the presentation matrix to account for deleted rows
    
    std::vector<index_t> rows;
    for(auto& entry : row_grade_map){
        rows.push_back(entry.first);
    }
    sort(rows.begin(), rows.end());
    hash_map<size_t, size_t> row_index_map;
    for(size_t i=0; i<rows.size(); i++){
        row_index_map[rows[i]] = i;
    }
    
    Matrix minimized_presentation;
    for(size_t i=0;i<presentation_matrix.size(); i++){
        SignatureColumn column(presentation_matrix[i].grade, presentation_matrix[i].signature_index);
        for(auto& entry : presentation_matrix[i]){
            if(row_grade_map.find(entry) != row_grade_map.end()){
                column.push(column_entry_t(entry, row_index_map[entry]));
            }
        }
        minimized_presentation.push_back(column);
    }
    
    std::cout << "Finished computing presentation of size "<< minimized_presentation.size();
    return std::pair<Matrix, hash_map<size_t, grade_t>>(minimized_presentation, row_grade_map);
}

std::pair<Matrix, hash_map<size_t, grade_t>> compute_presentation_schreyer(Matrix& image_generators, Matrix& kernel_generators, bool debug=false){
    /** Input a minimal set of generators for the image and a Groebner basis for the kernel.
     
     */
    std::cout << "Starting to compute presentation using Schreyer's algorithm"<< std::endl;
    
    sort(kernel_generators.begin(), kernel_generators.end(),[ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade < rhs.grade ;
         });
    
    for(size_t i=0; i<kernel_generators.size(); i++){
        kernel_generators[i].signature_index = i;
        kernel_generators[i].syzygy = SyzColumn();
        kernel_generators[i].syzygy.push(column_entry_t(1, i));
    }
    
    /* Translate image matrix */
    hash_map<index_t, std::vector<size_t>> pivots_indices;
    pivots_indices.reserve(kernel_generators.size());
    for(auto& column : kernel_generators){
        if(pivots_indices.find(column.get_pivot_index()) == pivots_indices.end()){
            pivots_indices[column.get_pivot_index()] = std::vector<size_t>();
        }
        pivots_indices[column.get_pivot_index()].push_back(column.signature_index);
    }
    
    Matrix translated_image_columns;
    for(auto& column : image_generators){
        SignatureColumn working_column(column);
        working_column.syzygy = SyzColumn();
        index_t pivot = working_column.get_pivot_index();
        while(pivot != -1){
            bool has_eliminated = false;
            std::vector<SignatureColumn> pivot_columns;
            for(auto& pivot_index : pivots_indices[pivot]){
                pivot_columns.push_back(kernel_generators[pivot_index]);
                if(kernel_generators[pivot_index].grade.leq_poset(working_column.grade)){
                    working_column.plus(kernel_generators[pivot_index]);
                    has_eliminated = true;
                    break;
                }
            }
            if(!has_eliminated){
                throw "Cannot express image column in terms of kernel.";
            }
            pivot = working_column.get_pivot_index();
        }
        translated_image_columns.push_back(SignatureColumn(working_column.grade, working_column.signature_index, working_column.syzygy));
    }
    
    Matrix syzygy_module = compute_syzygy_module(kernel_generators);
    std::cout << "Computed Syzygy module of size: "<< syzygy_module.size()<< std::endl;
    for(auto& column : translated_image_columns){
        syzygy_module.push_back(column);
    }
    
    std::pair<Matrix, Matrix> m_pair = compute_minimal_generating_set2(kernel_generators);
    Matrix& generating_set_kernel = m_pair.first;
    Matrix& change_of_basis_map = m_pair.second;
    
    std::cout << "Minimal generating set for kernel of size: "<< generating_set_kernel.size()<< std::endl;
    
    hash_map<size_t, grade_t> row_grade_map;
    for(size_t j=0; j<generating_set_kernel.size(); j++){
        row_grade_map[j] = generating_set_kernel[j].grade;
    }

    Matrix presentation;
    for(auto& syz_column : syzygy_module){
        SignatureColumn column(syz_column.grade, syz_column.signature_index);
        index_t pivot = syz_column.get_pivot().get_index();
        while(pivot != -1){
            column.plus(change_of_basis_map[pivot]);
            syz_column.pop_pivot();
            pivot = syz_column.get_pivot().get_index();
        }
        if(debug){
            if(syz_column.grade != column.grade){
                throw "Grade of columns has changed.";
            }
            SignatureColumn ker_column(syz_column.grade, syz_column.signature_index);
            SignatureColumn copy(column);
            pivot = copy.get_pivot().get_index();
            while(pivot != -1){
                ker_column.plus(generating_set_kernel[pivot]);
                copy.pop_pivot();
                pivot = copy.get_pivot().get_index();
            }
            if(ker_column.get_pivot().get_index() != -1){
                copy.print();
                std::cout << "\n";
                column.print();
                std::cout << "\n";
                ker_column.print();
                throw "Column non-zero";
            }
        }
        if(column.get_pivot().get_index() != -1){
            presentation.push_back(column);
        }
    }

    std::cout << "Computed non-minimal presentation of size: " << presentation.size()<< std::endl;
    
    // Compute minimal generating set from syz-module and image columns
    sort(presentation.begin(), presentation.end(),[ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade < rhs.grade ;
         });
    
    hash_map<size_t, size_t> pairs = compute_local_pairs(presentation, row_grade_map);
    
    std::set<size_t> negative_columns;
    presentation = compute_global_columns(presentation, pairs, negative_columns);
    
    for(auto& entry : pairs){
        if(row_grade_map.find(entry.first) != row_grade_map.end()){
            row_grade_map.erase(entry.first);
        }
    }
    
    std::vector<index_t> rows;
    for(auto& entry : row_grade_map){
        rows.push_back(entry.first);
    }
    sort(rows.begin(), rows.end());
    hash_map<size_t, size_t> row_index_map;
    for(size_t i=0; i<rows.size(); i++){
        row_index_map[rows[i]] = i;
    }
    
    for(size_t i=0;i<presentation.size(); i++){
        SignatureColumn column(presentation[i].grade, presentation[i].signature_index);
        for(auto& entry : presentation[i]){
            if(row_grade_map.find(entry) != row_grade_map.end()){
                column.push(column_entry_t(entry, row_index_map[entry]));
            }
        }
        presentation[i].swap(column);
    }
    presentation = compute_minimal_generating_set(presentation);
    return std::pair<Matrix, hash_map<size_t, grade_t>>(presentation, row_grade_map);
}


/* Python interface */

struct GradedMatrix{
    std::vector< std::vector<std::pair<int, int>> > matrix;
    std::vector< std::vector<int> > column_grades;
    std::vector< std::vector<int> > row_grades;
};

Metric* parse_metric(int metric_index){
    switch(metric_index){
        case 0:
            return new SquaredEuclideanMetric();
        default:
            throw "Failed to parse metric";
    }
}

Filter* parse_filter(int filter_index){
    if(filter_index >= 0){
        return new XFilter(filter_index);
    }
    switch(filter_index){
        case -1:
            return new XFilter(0); //TODO: implement more filters.
        default:
            throw "Failed to parse metric";
    }
}

Matrix translateInputMatrix(std::vector<std::vector<int>>& matrix, std::vector<std::vector<int>>& column_grades){
    Matrix M;
    for(size_t i=0; i<matrix.size(); i++){
        grade_t grade;
        for(size_t j=0; j<column_grades[i].size(); j++){
            grade.push_back(column_grades[i][j]);
        }
        SignatureColumn column(grade, i);
        for(size_t j=0; j<matrix[i].size(); j++){
            if(matrix[i][j] != 0){
                column.push(column_entry_t(1, j)); // TODO: implement support for other modulus than Z/2Z.
            }
        }
        column.syzygy.push(column_entry_t(1, i));
        M.push_back(column);
    }
    return M;
}

GradedMatrix translateOutputMatrix(Matrix& matrix){
    GradedMatrix graded_matrix;
    for(SignatureColumn column : matrix){
        std::vector<std::pair<int, int>> sparse_column;
        while(!column.empty()){
            column_entry_t entry = column.pop_pivot();
            if(entry.get_index() != -1){
                sparse_column.push_back(std::pair<int, int>(entry.get_value(), entry.get_index()));
            }
        }
        graded_matrix.matrix.push_back(sparse_column);
        std::vector<int> column_grade;
        for(size_t grade_index=0; grade_index<column.grade.size(); grade_index++){
            column_grade.push_back((int)column.grade[grade_index]);
        }
        graded_matrix.column_grades.push_back(column_grade);
    }
    return graded_matrix;
}

GradedMatrix compute_kernel(std::vector<std::vector<int>>& matrix, std::vector<std::vector<int>>& row_grades, std::vector<std::vector<int>>& column_grades){
    Matrix M = translateInputMatrix(matrix, column_grades);
    Matrix kernel;
    if(column_grades.size()>0 && column_grades[0].size()==2){
        kernel = computeKernel_2p(M);
    }else{
        std::pair<Matrix, Matrix> gbs = computeGroebnerBases_gradeopt_min(M);
        kernel = gbs.first;
    }
    GradedMatrix Ker = translateOutputMatrix(kernel);
    Ker.row_grades = column_grades;
    return Ker;
}

std::pair<GradedMatrix, GradedMatrix> groebner_bases(std::vector<std::vector<int>>& matrix, std::vector<std::vector<int>>& row_grades, std::vector<std::vector<int>>& column_grades){
    Matrix M = translateInputMatrix(matrix, column_grades);
    std::pair<Matrix, Matrix> gbs = computeGroebnerBases_gradeopt_min(M);
    GradedMatrix Im = translateOutputMatrix(gbs.first);
    GradedMatrix Ker = translateOutputMatrix(gbs.second);
    Ker.row_grades = column_grades;
    Im.row_grades = row_grades;
    return std::pair<GradedMatrix, GradedMatrix>(Im, Ker);
}

GradedMatrix presentation_FIrep(std::vector<std::vector<int>>& high_matrix, std::vector<std::vector<int>>& column_grades_h, std::vector<std::vector<int>>& low_matrix, std::vector<std::vector<int>>& column_grades_l){
    Matrix M_h = translateInputMatrix(high_matrix, column_grades_h);
    Matrix M_l = translateInputMatrix(low_matrix, column_grades_l);
    auto start = std::chrono::high_resolution_clock::now();
    std::pair<Matrix, hash_map<size_t, grade_t>> presentation_output;
    if(M_l.size()>0 && M_l[0].grade.size()==2){
        Matrix kernel = computeKernel_2p(M_l);
        presentation_output = compute_presentation_2p(M_h, kernel);
    }else{
        presentation_output = computePresentationDeg_imopt(M_h, M_l);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    
    std::cout << "Time elapsed: " << duration.count() << " seconds" << std::endl;
    GradedMatrix graded_matrix = translateOutputMatrix(presentation_output.first);
    for(std::pair<size_t, grade_t> entry : presentation_output.second){
        std::vector<int> row_grade;
        for(size_t grade_index=0; grade_index<entry.second.size(); grade_index++){
            row_grade.push_back((int)entry.second[grade_index]);
        }
        graded_matrix.row_grades.push_back(row_grade);
    }
    return graded_matrix;
}

GradedMatrix presentation_dm(std::vector<std::vector<std::vector<input_t>>>& distance_matrices, std::vector<input_t>& max_metric_values, std::vector<std::vector<input_t>>& filters, int hom_dim){
    std::pair<Matrix, Matrix> boundary_matrices = compute_boundary_matrices_dm(distance_matrices, max_metric_values, filters, hom_dim);
    verify_kernel(boundary_matrices.first, boundary_matrices.second);
    for(size_t i=0; i<boundary_matrices.second.size(); i++){
        boundary_matrices.second[i].syzygy.push(column_entry_t(1, i));
    }
    Matrix input_copy;
    for(auto& column : boundary_matrices.second){
        input_copy.push_back(SignatureColumn(column.grade, 0, column));
    }
    std::pair<Matrix, hash_map<size_t, grade_t>> presentation;
    if(boundary_matrices.second.size()>0 && boundary_matrices.second[0].grade.size()==2){
        Matrix kernel = computeKernel_2p(boundary_matrices.second);
        presentation = compute_presentation_2p(boundary_matrices.first, kernel);
    }else{
        std::pair<Matrix, Matrix> gbs = computeGroebnerBases(boundary_matrices.second);
        presentation = computePresentationDeg_imopt(boundary_matrices.first, gbs.second);
    }
    GradedMatrix graded_matrix = translateOutputMatrix(presentation.first);
    for(std::pair<size_t, grade_t> entry : presentation.second){
        std::vector<int> row_grade;
        for(size_t grade_index=0; grade_index<entry.second.size(); grade_index++){
            row_grade.push_back((int)entry.second[grade_index]);
        }
        graded_matrix.row_grades.push_back(row_grade);
    }
    return graded_matrix;
}

GradedMatrix presentation(std::vector<std::vector<input_t>>& _points, std::vector<int>& _metrics, std::vector<input_t>& _max_metric_values, std::vector<int>& _filters, int hom_dim){
    std::vector<Metric*> metrics;
    for(size_t i=0; i<_metrics.size(); i++){
        metrics.push_back(parse_metric(_metrics[i]));
    }
    std::vector<Filter*> filters;
    for(size_t i=0; i<filters.size(); i++){
        filters.push_back(parse_filter(_filters[i]));
    }
    std::vector<input_t> max_metric_values;
    for(size_t i=0; i<_metrics.size(); i++){
        max_metric_values.push_back(_max_metric_values[i]);
    }
    std::pair<Matrix, Matrix> boundary_matrices = compute_boundary_matrices(_points, metrics, filters, max_metric_values, hom_dim);
    verify_kernel(boundary_matrices.first, boundary_matrices.second);
    for(auto& metric : metrics){
        delete metric;
    }
    for(auto& filter : filters){
        delete filter;
    }
    for(size_t i=0; i<boundary_matrices.second.size(); i++){
        boundary_matrices.second[i].syzygy.push(column_entry_t(1, i));
    }
    Matrix input_copy;
    for(auto& column : boundary_matrices.second){
        input_copy.push_back(SignatureColumn(column.grade, 0, column));
    }
    std::pair<Matrix, hash_map<size_t, grade_t>> presentation;
    if(boundary_matrices.second.size()>0 && boundary_matrices.second[0].grade.size()==2){
        Matrix kernel = computeKernel_2p(boundary_matrices.second);
        presentation = compute_presentation_2p(boundary_matrices.first, kernel);
    }else{
        presentation = computePresentationDeg_imopt(boundary_matrices.first, boundary_matrices.second);
    }
    GradedMatrix graded_matrix = translateOutputMatrix(presentation.first);
    for(std::pair<size_t, grade_t> entry : presentation.second){
        std::vector<int> row_grade;
        for(size_t grade_index=0; grade_index<entry.second.size(); grade_index++){
            row_grade.push_back((int)entry.second[grade_index]);
        }
        graded_matrix.row_grades.push_back(row_grade);
    }
    return graded_matrix;
}


/*
 Testing
 */

std::vector<double> run_gb_singatures_pres(Matrix& input_matrix, Matrix& image, bool debug=false){
    auto start = std::chrono::high_resolution_clock::now();
    std::pair<Matrix, hash_map<size_t, grade_t>> presentation_sign = computePresentationDeg_imopt(image, input_matrix, debug);
    auto stop = std::chrono::high_resolution_clock::now();
    auto b_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << " Presentation size: " << presentation_sign.first.size() << std::endl;//gbs2.size() << std::endl;
    double time = b_time.count();
    double size = presentation_sign.first.size();
    return std::vector<double>{time, res_memory, size};
}

std::vector<double> run_gb_shreyer_pres(Matrix& input_matrix, Matrix& image, bool debug=false){
    auto start = std::chrono::high_resolution_clock::now();
    Matrix kernel = computeGroebnerBases_gradeopt(input_matrix).second;
    std::pair<Matrix, hash_map<size_t, grade_t>> presentation = compute_presentation_schreyer(image, kernel);
    auto stop = std::chrono::high_resolution_clock::now();
    auto b_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << " Presentation size: " << presentation.first.size() << std::endl;//gbs2.size() << std::endl;
    double time = b_time.count();
    double size = presentation.first.size();
    return std::vector<double>{time, res_memory, size};
}

void critical_points_geometric_pres(){
    get_mem_usage(virt_memory, res_memory);
    double start_mem = res_memory;
    double metric_thresh = 3.5;
    int n = 10;
    std::vector<std::vector<double>> times;
    std::vector<std::vector<double>> resident_memory;
    std::vector<std::vector<double>> output_size;
    std::vector<double> boundary_matrix_size;
    while(n < 201){
        std::cout << "\n Iteration: " << n << std::endl;
        uint32_t seed = 1;
        std::vector<std::vector<input_t>> points = time_varying_point_cloud(10, 10*n, 3, seed);//time_varying_point_cloud(60, 100, 3, 1);
        std::vector<std::vector<input_t>> xpoints = points;
        for(auto& p:xpoints){
            p.pop_back();
        }
        
        std::vector<Metric*> metrics;
        metrics.push_back(new SquaredEuclideanMetric());
        std::vector<Filter*> filters;
        filters.push_back(new XFilter(3, 1));
        filters.push_back(new XFilter(3, -1));
        //filters.push_back(new XFilter(0, 1));
        //filters.push_back(new XFilter(0, -1));
        //filters.push_back(new XFilter(1, 1));
        //filters.push_back(new XFilter(1, -1));
        //filters.push_back(new XFilter(2, 1));
        //filters.push_back(new XFilter(2, -1));
        //filters.push_back(new DensityFilter(xpoints));
        
        std::vector<input_t> max_metric_values;
        max_metric_values.push_back(metric_thresh);
        
        std::pair<Matrix, Matrix> boundaries = compute_boundary_matrices(points, metrics, filters, max_metric_values, 1);
        
        for(size_t i=0; i<boundaries.second.size(); i++){
            boundaries.second[i].signature_index = i;
        }
        
        sort(boundaries.second.begin(), boundaries.second.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
             {
                 return lhs.grade.lt_colex(rhs.grade);
             });
        
        //Reindex image columns
        hash_map<size_t, size_t> map;
        for(size_t i=0; i<boundaries.second.size();i++){
            map[boundaries.second[i].signature_index] = i;
        }
        
        for(size_t i=0; i<boundaries.first.size(); i++){
            SyzColumn c;
            while(boundaries.first[i].get_pivot().get_index() != -1){
                c.push(column_entry_t(1, map[boundaries.first[i].get_pivot().get_index()]));
                boundaries.first[i].pop_pivot();
            }
            while(c.get_pivot().get_index() != -1){
                boundaries.first[i].push(c.get_pivot());
                c.pop_pivot();
            }
        }
        
        sort(boundaries.first.begin(), boundaries.first.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
             {
                 return lhs.grade < rhs.grade;
             });
        
        Matrix boundary = boundaries.second;
        Matrix image = compute_minimal_generating_set(boundaries.first);
        
        std::cout << "Size of boundary matrix: " << boundary.size() << " Size of image matrix: " << image.size() << std::endl;
        
        for(size_t i=0; i<boundary.size(); i++){
            boundary[i].signature_index = i;
            boundary[i].syzygy = SyzColumn();
            boundary[i].syzygy.push(column_entry_t(1, i));
        }
        
        //Test is image maps into kernel
        for(auto& column : image){
            SignatureColumn c(column);
            SignatureColumn z(column.grade, -1);
            while(c.get_pivot_index() != -1){
                z.plus(boundary[c.get_pivot_index()]);
                c.pop_pivot();
            }
            if(z.get_pivot_index() != -1){
                throw "Image does not map to kernel";
            }
        }
        
        Matrix input_copy;
        for(auto& c : boundary){
            input_copy.push_back(SignatureColumn(c));
        }
        Matrix image_copy;
        for(auto& c : image){
            image_copy.push_back(SignatureColumn(c));
        }
        std::vector<double> gb_res = run_gb_singatures_pres(input_copy, image_copy);
        input_copy.clear();
        for(auto& c : boundary){
            input_copy.push_back(SignatureColumn(c));
        }
        image_copy.clear();
        for(auto& c : image){
            image_copy.push_back(SignatureColumn(c));
        }
        std::vector<double> gb_schreyer = run_gb_shreyer_pres(input_copy, image_copy);
        
        times.push_back(std::vector<double>{gb_schreyer[0], gb_res[0]});
        resident_memory.push_back(std::vector<double>{gb_schreyer[1]-start_mem, gb_res[1]-start_mem});
        output_size.push_back(std::vector<double>{gb_schreyer[2], gb_res[2]});
        boundary_matrix_size.push_back(boundary.size());
        
        std::cout << "$"<< metric_thresh << "$ & " << times[times.size()-1][0]/1000 << " & " << resident_memory[times.size()-1][0] << " & " << times[times.size()-1][1]/1000 << " & " << resident_memory[times.size()-1][1]  << "\n";
        
        n+= 10;
    }
    for(size_t i=0; i<times.size(); i++){
        std::cout << "$"<< metric_thresh+i*0.05 << "$ & " << times[i][0]/1000 << " & " << resident_memory[i][0] << " & " << times[i][1]/1000 << " & " << resident_memory[i][1] << "\n";
    }
    std::cout << "\n";
    std::cout << "\n Times GBS+Schreyer \n";
    for(auto& v : times){
        std::cout << v[0] << ", ";
    }
    std::cout << "\n Times PresentationPair \n";
    for(auto& v : times){
        std::cout << v[1] << ", ";
    }
    std::cout << "\n Resident memory GBS+Schreyer\n";
    for(auto& v : resident_memory){
        std::cout << v[0] << ", ";
    }
    std::cout << "\n Resident memory PresentationPair\n";
    for(auto& v : resident_memory){
        std::cout << v[1] << ", ";
    }
    
    std::cout << "\n Input size \n";
    for(auto& v : boundary_matrix_size){
        std::cout << v << ", ";
    }
}


/* Input */

void print_usage_and_exit(int exit_code) {
    std::cerr
    << "Usage: "
    << "mph "
    << "[options] [filename]" << std::endl
    << std::endl
    << "Options:" << std::endl
    << std::endl
    << "  --help           print this screen" << std::endl
    << "  --dim <k>        compute presentation matrix of the k-th persistent homology module" << std::endl
    << "  --firep <k>      input file with rivet <firep> file format" << std::endl
    << std::endl;
    exit(exit_code);
}


int main(int argc, char** argv) {
    const char* filename = nullptr;
    
    int dim_max = 1;
    bool firep = false;
    
    if(argc==1){
        print_usage_and_exit(0);
    }
    
    for (index_t i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--help") {
            print_usage_and_exit(0);
        } else if (arg == "--dim") {
            std::string parameter = std::string(argv[++i]);
            size_t next_pos;
            dim_max = (int)std::stol(parameter, &next_pos);
            if (next_pos != parameter.size()) print_usage_and_exit(-1);
        } else if (arg == "--firep"){
            firep = true;
        } else {
            if (filename) { print_usage_and_exit(-1); }
            filename = argv[i];
        }
    }
    
    std::pair<Matrix, hash_map<size_t, grade_t>> presentation;
    if(firep){
        Matrix high_matrix, low_matrix;
        std::ifstream file_stream(filename);
        if (file_stream.fail()) {
            std::cerr << "couldn't open file " << filename << std::endl;
            exit(-1);
        }
        read_input_file<SyzColumn, SyzColumn>(file_stream, high_matrix, low_matrix);
        presentation = computePresentationDeg_imopt(high_matrix, low_matrix);
    } else{
        std::ifstream file_stream(filename);
        if (filename && file_stream.fail()) {
            std::cerr << "couldn't open file " << filename << std::endl;
            exit(-1);
        }
        
        // TODO: how should choice of metrics and filters be specified in input?
        
        std::vector<Metric*> metrics;
        metrics.push_back(new SquaredEuclideanMetric());
        std::vector<Filter*> filters;
        filters.push_back(new Filter());
        std::vector<input_t> max_metric_values;
        max_metric_values.push_back(10);
        std::vector<std::vector<input_t>> points = read_point_cloud(file_stream, 0);
        std::pair<Matrix, Matrix> boundary_matrices = compute_boundary_matrices(points, metrics, filters, max_metric_values, dim_max);
        verify_kernel(boundary_matrices.first, boundary_matrices.second);
        for(auto& metric : metrics){
            delete metric;
        }
        for(auto& filter : filters){
            delete filter;
        }
        for(size_t i=0; i<boundary_matrices.second.size(); i++){
            boundary_matrices.second[i].syzygy.push(column_entry_t(1, i));
        }
        presentation = computePresentationDeg_imopt(boundary_matrices.first, boundary_matrices.second);
    }
    
    std::cout << "Presentation matrix: \n";
    presentation.first.print();
    
    std::cout << "Column grades: \n";
    for(auto& column: presentation.first){
        column.grade.print();
    }
    std::cout << "Row grades: \n";
    for(auto& it: presentation.second){
        it.second.print();
    }
    exit(0);
}
