#=================================#
# This file contains all the
# "control" functions
#=================================#

# Possible elements of myAlerts:
#

myCallAlerts = function(call, name, myType, nParents=1, mustBeThere=FALSE, prefix=""){
    # This function basically calls the function myAlerts, but the arguments are different

    if( name %in% names(call) ){
        # we check the element exists => to provide a fine error
        what = call[[name]]
        val = try(eval.parent(what, nParents), silent = TRUE)
        # browser()
        if(inherits(val, "try-error") ){
            if(inherits(what, "name")){
                # it means the variable was not found
                stop(prefix,"For argument '",name,"': object '",what,"' not found.", call. = FALSE)
            } else {
                stop(prefix,"For argument '",name,"': expression ",as.character(as.expression(what))," could not be evaluated.", call. = FALSE)
            }

        } else {
            a = myAlerts(val, name, myType, prefix)
            return(a)
        }
    } else if(mustBeThere) {
        stop(prefix, "The argument '", name, "' must be provided.", call. = FALSE)
    }
}

myAlerts = function(x, name, myType, prefix="", charVec){
    # Format of my types:
    #   - single => must be of lenght one
    #   - Vector => must be a vector
    #   - Matrix => must be a matrix
    #   - GE/GT/LE/LT: greater/lower than a given value
    #   - predefinedType => eg: numeric, integer, etc
    #   - match.arg => very specific => should match the charVec
    # If there is a parenthesis => the class must be of specified types:
    # ex: "(list, data.frame)" must be a list of a data.frame
	
	ignore.case = TRUE

    firstMsg = paste0(prefix,"The argument '",name,"' ")

    # simple function to extract a pattern
    # ex: if my type is VectorIntegerGE1 => myExtract("GE[[:digit:]]+","VectorIntegerGE1") => 1
    myExtract = function(expr, text, trim=2){
        start = gregexpr(expr,text)[[1]] + trim
        length = attr(start,"match.length") - trim
        res = substr(text,start,start+length-1)
        as.numeric(res)
    }

    #
    # General types handling
    #

    loType = tolower(myType)

    if(grepl("single",loType)){
        if(length(x)!=1) stop(firstMsg,"must be of length one.", call. = FALSE)
    }

    if(grepl("vector",loType) && !grepl("factor",loType)){
        if(!is.vector(x)) stop(firstMsg,"must be a vector.", call. = FALSE)
        if(is.list(x)) stop(firstMsg,"must be a vector (and not a list).", call. = FALSE)
    }
    
    res = checkTheTypes(loType, x)
    if(!res$OK) stop(firstMsg,res$message, call. = FALSE)

    # # INTEGER is a restrictive type that deserves some explanations (not included in getTheTypes)
    # if(grepl("integer",loType)){
    #     if(grepl("single",loType)){
    #         if(!is.numeric(x)) stop(firstMsg,"must be an integer (right now it is not even numeric).", call. = FALSE)
    #         if(!(is.integer(x) || x%%1==0)) stop(firstMsg,"must be an integer.", call. = FALSE)
    #     } else {
    #         if(!is.numeric(x)) stop(firstMsg,"must be composed of integers (right now it is not even numeric).", call. = FALSE)
    #         if(!(is.integer(x) || all(x%%1==0))) stop(firstMsg,"must be composed of integers.", call. = FALSE)
    #     }
    # }

    # GE: greater or equal // GT: greater than // LE: lower or equal // LT: lower than
    if(grepl("ge[[:digit:]]+",loType)){
        n = myExtract("ge[[:digit:]]+", loType)
        if( !all(x>=n) ) stop(firstMsg,"must be greater than, or equal to, ", n,".", call. = FALSE)
    }
    if(grepl("gt[[:digit:]]+",loType)){
        n = myExtract("gt[[:digit:]]+", loType)
        if( !all(x>n) ) stop(firstMsg,"must be strictly greater than ", n,".", call. = FALSE)
    }
    if(grepl("le[[:digit:]]+",loType)){
        n = myExtract("le[[:digit:]]+", loType)
        if( !all(x<=n) ) stop(firstMsg,"must be lower than, or equal to, ",n,".", call. = FALSE)
    }
    if(grepl("lt[[:digit:]]+",loType)){
        n = myExtract("lt[[:digit:]]+", loType)
        if( !all(x<n) ) stop(firstMsg,"must be strictly lower than ", n,".", call. = FALSE)
    }

    #
    # Specific Types Handling
    #

    if(grepl("match.arg",loType)){
    	if(ignore.case){
    		x = toupper(x)
    		newCharVec = toupper(charVec)
    	} else {
    		newCharVec = charVec
    	}
    	
        if( is.na(pmatch(x, newCharVec)) ){
            n = length(charVec)
            if(n == 1){
                msg = paste0("'",charVec,"'")
            } else {
                msg = paste0("'", paste0(charVec[1:(n-1)], collapse="', '"), "' or '",charVec[n],"'")
            }
            stop(firstMsg, "must be one of:\n", msg, ".", call. = FALSE)
        } else {
        	qui = pmatch(x, newCharVec)
        	return(charVec[qui])
        }
    }
}

matchTypeAndSetDefault = function(myList, myDefault, myTypes, prefix){
    # Cette fonction:
    #   i) verifie que tous les elements de la liste sont valides
    #   ii) mes les valeurs par defauts si elles certaines valeurs sont manquantes
    #   iii) Envoie des messages d'erreur si les typages ne sont pas bons
    # En fait cette fonction "coerce" myList en ce qu'il faudrait etre (donne par myDefault)

    # 1) check that the names of the list are valid
    if(is.null(myList)) myList = list()
    list_names = names(myList)

    if(length(list_names)!=length(myList) || any(list_names=="")){
        stop(prefix,"The elements of the list should be named.", call. = FALSE)
    }

    obj_names = names(myDefault)

    isHere = pmatch(list_names,obj_names)

    if(anyNA(isHere)){
        if(sum(is.na(isHere))==1) stop(prefix, "The following argument is not defined: ",paste(list_names[is.na(isHere)],sep=", "), call. = FALSE)
        else stop(prefix, "The following arguments are not defined: ",paste(list_names[is.na(isHere)],sep=", "), call. = FALSE)
    }

    # 2) We set the default values and run Warnings
    res = list()
    for(i in 1:length(obj_names)){
        obj = obj_names[i]
        qui = which(isHere==i) # qui vaut le numero de l'objet dans myList
        type = myTypes[i] # we extract the type => to control for "match.arg" type
        if(length(qui)==0){
            # we set to the default if it's missing
            if(type == "match.arg") {
                res[[obj]] = myDefault[[i]][1]
            } else {
                res[[obj]] = myDefault[[i]]
            }
        } else {
            # we check the integrity of the value
            val = myList[[qui]]
            if(type == "match.arg"){
                # If the value is to be a match.arg => we use our controls and not
                # directly the one of the function match.arg()
                charVec = myDefault[[i]]
                myAlerts(val, obj, "singleCharacterMatch.arg", prefix, charVec)
                val = match.arg(val, charVec)
            } else {
                myAlerts(val, obj, type, prefix)
            }

            res[[obj]] = val
        }
    }

    return(res)
}



checkTheTypes = function(str, x){
	# This function takes in a character string describing the types of the 
	# element x => it can be of several types 
	
	# types that are controlled for:
	allTypes = c("numeric", "integer", "character", "logical", "list", "data.frame", "matrix", "factor")
	
	OK = FALSE
	message = c()
	
	for(type in allTypes){
		
		if(grepl(type, str)){
			
			# we add the type of the control
			message = c(message, type)
			
			if(type == "numeric"){
				if(!OK & is.numeric(x)){
					OK = TRUE
				}
			} else if(type == "integer"){
				if(is.numeric(x) && (is.integer(x) || all(x%%1==0))){
					OK = TRUE
				}
			} else if(type == "character"){
				if(is.character(x)){
					OK = TRUE
				}
			} else if(type == "logical"){
				if(is.logical(x)){
					OK = TRUE
				}
			} else if(type == "list"){
				if(is.list(x)){
					OK = TRUE
				}
			} else if(type == "data.frame"){
				if(is.data.frame(x)){
					OK=TRUE
				}
			} else if(type == "matrix"){
				if(is.matrix(x)){
					OK = TRUE
				}
			} else if(type == "factor"){
				if(is.factor(x)){
					OK = TRUE
				}
			}  
		}
		
		if(OK) break
	}
	
	if(length(message) == 0) OK = TRUE #ie there is no type to be searched
	else if(length(message) >= 3){
		n = length(message)
		message = paste0("must be of type: ",  paste0(message[1:(n-1)], collapse = ", "), " or ", message[n], ".")
	} else {
		message = paste0("must be of type: ",  paste0(message, collapse = " or "), ".")
	}
	
	
	return(list(OK=OK, message=message))
}

