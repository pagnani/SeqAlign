module SeqAlign

export smithwaterman, needlemanwunsch

#                  A    C   C   T
# const ScoreDNA = [[0.6 0.1 0.2 0.1]; # A 
#                   [0.1 0.6 0.1 0.2]; # C
#                   [0.2 0.1 0.6 0.1]; # G
#                   [0.1 0.2 0.1 0.6]; # T
#                   ]
const ScoreDNA = [[ 5  -3 -3 -3]; # A 
                  [-3   5 -3 -3]; # C
                  [-3  -3  5 -3]; # G
                  [-3  -3 -3  5]; # T
                  ]

const GapScore = -4

# const ScoreDNA = [[ 1  -1 -1 -1]; # A 
#                   [-1   1 -1 -1]; # C
#                   [-1  -1  1 -1]; # G
#                   [-1  -1 -1  1]; # T
#                   ]


function needlemanwunsch(seq1::String, seq2::String)

    IterMat = CreateMatrixNW(seq1,seq2)
    align =  BackTrackNW( IterMat, seq1, seq2)    
    stringal = SimplePrintAlign(align)
    
    return stringal 
    
end

function smithwaterman(seq1::String, seq2::String)
    
    IterMat = CreateMatrixSW(seq1,seq2)
    align =  BackTrackSW( IterMat, seq1, seq2)    
    stringal = SimplePrintAlign(align)
    
     return stringal
end

function CreateMatrixNW(seq1::String,seq2::String)

    N1 = length(seq1)
    N2 = length(seq2)
    iseq1 = TransformSeq(seq1)
    iseq2 = TransformSeq(seq2)


    IterMat = zeros(eltype(ScoreDNA),N1+1, N2+1)
    @inbounds begin
        for i=1:N2+1
            IterMat[1,i] = -(i-1)*GapScore
        end
        for i=1:N1+1
            IterMat[i,1] = -(i-1)*GapScore
        end
    end
    loopmat!(IterMat,iseq1,iseq2, N1,N2)
    return  IterMat
end



function loopmat!(IterMat::Array{eltype(ScoreDNA),2},iseq1::Array{Int,1},iseq2::Array{Int,1}, N1::Int, N2::Int)
    @inbounds begin
        for j = 2:N2+1
           for i=2:N1+1
               scra1 = IterMat[i-1,j-1] + ScoreDNA[iseq1[i-1],iseq2[j-1]]
               scra2 = IterMat[i-1,j  ] + GapScore
               scra3 = IterMat[i  ,j-1] + GapScore 
               IterMat[i,j] = max(scra1,scra2,scra3)
           end
        end
    end
end



function BackTrackNW( IterMat::Array{eltype(ScoreDNA),2}, seq1::String, seq2::String)

    siz = size(IterMat)
    iseq1 = TransformSeq(seq1)
    iseq2 = TransformSeq(seq2)

    lalign  = (Char,Char,Char)[]
    align = Array(typeof(lalign),1)
    i = siz[1]
    j = siz[2]
    
    while i > 1 && j > 1
        if i > 1 && j > 1 && (IterMat[i,j] == IterMat[i-1,j-1] + ScoreDNA[iseq1[i-1],iseq2[j-1]])
            if seq1[i-1] == seq2[j-1]
                push!(lalign, (seq1[i-1],seq2[j-1],'+'))
            else
                push!(lalign, (seq1[i-1],seq2[j-1],'o'))
            end
            i -= 1
            j -= 1
        elseif i > 1 && (IterMat[i,j] == IterMat[i-1,j] + GapScore)
            push!(lalign,(seq1[i-1],'-','-'))
            i -= 1
        elseif j > 1 && (IterMat[i,j] == IterMat[i,j-1] + GapScore)
            push!(lalign,('-',seq2[j-1],'-'))
            j -= 1
        else
            error("What am I doing here")
        end     
    end 
    
    for l=i-1:-1:1
        push!(lalign,(seq1[l],'-','-'))
    end
    for l=j-1:-1:1
        push!(lalign,('-',seq2[l],'-'))
    end

    align[1] = copy(lalign)
    return align
end


function BackTrackSW( IterMat::Array{eltype(ScoreDNA),2}, seq1::String, seq2::String)

    iseq1 = TransformSeq(seq1)
    iseq2 = TransformSeq(seq2)
    

    siz = size(IterMat)
    totalmassimo = maximum(IterMat)
    pairmax = find(IterMat .== totalmassimo)
    nummatch = length(pairmax)
    
    lalign  = (Char,Char,Char)[]
    align = Array(typeof(lalign),nummatch)



    for t = 1:nummatch
        i, j  = ind2sub(siz, pairmax[t])
        while i > 1 && j > 1 && IterMat[i,j] > 0 
            if i > 1 && j > 1 && ( IterMat[i,j] == IterMat[i-1,j-1] + ScoreDNA[iseq1[i-1],iseq2[j-1]])
                if seq1[i-1] == seq2[j-1]
                    push!(lalign, (seq1[i-1],seq2[j-1],'+'))
                else
                    push!(lalign, (seq1[i-1],seq2[j-1],'o'))
                end           
                i -= 1
                j -= 1
            elseif  i > 1 && (IterMat[i,j] == IterMat[i-1,j] + GapScore) 
                push!(lalign,(seq1[i-1],'-','-'))
                i -= 1
            elseif j > 1 && (IterMat[i,j] == IterMat[i,j-1] + GapScore)
                push!(lalign,('-',seq2[j-1],'-'))
                j -= 1 
            else
                error("What am I doing here?")
            end
        end
        align[t] = copy(lalign)
        empty!(lalign)
    end
    return align
end


function CreateMatrixSW(seq1::String,seq2::String)

    iseq1 = TransformSeq(seq1)
    iseq2 = TransformSeq(seq2)


    N1 = length(iseq1)
    N2 = length(iseq2)


    IterMat = zeros(eltype(ScoreDNA),N1+1, N2+1)

    loopmat!(IterMat,iseq1, iseq2, N1,N2)

    return IterMat
end

function SimplePrintAlign(align::Array{Array{(Char,Char,Char),1},1})

    numal = length(align)

    
    stringal = Array((String, String, String),numal);

    for i=1:numal
        str1 = "";
        str2 = "";
        str3 = "";
        nal = length(align[i])
        for j = nal:-1:1
            str3 *= @sprintf("%c", align[i][j][3])
            str2 *= @sprintf("%c", align[i][j][2])
            str1 *= @sprintf("%c", align[i][j][1])
        end

        stringal[i] = (str1,str2,str3);
    end
    return stringal
end


function TransformSeq(seq::String)
    N = length(seq)
    vec = [ letter2num(seq[i]) for i=1:N ]
    return vec
end


let alphabetDNA = [1, 2, 3, 4, 5]
#                  A  C  G  T  ? 
    global letter2num
    function letter2num(c::Char)
        if c == 'A' || c == 'a'
            return 1
        elseif c == 'C' || c == 'c'
            return 2
        elseif c == 'G' || c == 'g'
            return 3
        elseif c == 'T' || c == 't' || c == 'U' || c == 'u'
            return 4
        else
            error("unknown base $c")
        end
    end
end   
end
