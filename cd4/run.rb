#!/usr/bin/env ruby
require 'optparse'
require 'fileutils'
require '/home/cew54/slurmer/Qsub.rb'
require 'pp'

## empty array for commands

OPTIONS = {}
OPTIONS[:int] = false
OPTIONS[:autorun] = false
OPTIONS[:all] = false
OPTIONS[:val] = true
OPTIONS[:pchic] = true
OptionParser.new do |opts|
  opts.banner = "Usage: runguess.rb [OPTIONS] COMMAND [DIR]"

  # opts.on("-a", "--[no-]all", "(Re-)run all outputs regardless of whether it already exists") do |i|
  #   OPTIONS[:all] = i
  # end
  
  opts.on("-i", "--[no-]interactive", "Run each job in serial on the interactive log in node") do |i|
    OPTIONS[:int] = i
  end
  
  opts.on("-r", "--autoRun", "Run each job on the queue without user input") do |i|
    OPTIONS[:autorun] = i
  end
  
  opts.on("-v", "--[no-]verbose", "Run verbosely") do |v|
    OPTIONS[:verbose] = v
  end

  opts.on("-h", "--help", "Show help") do |h|
    OPTIONS[:help] = h
  end

  opts.on("--pchic", "Just do pchic sets") do |v|
    OPTIONS[:val] = false
  end
  opts.on("--val", "Just do val sets") do |v|
    OPTIONS[:pchic] = false
  end
  # opts.on("-q", "--[no-]queue", "Use qCom.sh in output") do |v|
  #   OPTIONS[:queue] = v
  # end
  # opts.on("-n", "--nohup", "Prepend with nohup") do |v|
  #   OPTIONS[:nohup] = v
  # end
  # opts.on("-b", "--[no-]background", "Append with &") do |v|
  #   OPTIONS[:background] = v
  # end
end.parse!
COMMAND = ARGV.shift

def usage()
  puts OPTIONS
  puts "Usage: run.rb [OPTIONS] COMMAND

  COMMANDS are:
      brun  :: run reps in a block using --array
      run   :: run additional reps as needed
      summ  :: summarise current state
      corr  :: recalculate correlation
      mppc  :: store mppc for baits which pass the mppc >= 0.75
      clean :: remove slurm files      
  "
end
if OPTIONS[:help] then
  usage()
  exit 0
end
commands = []
COMFILE="tmprun.sh"
args = {:job=>COMMAND,
        :time=>"12:00:00",
        :tasks=>'1',
        :cpus=>'1',
	:account=>'CWALLACE-SL2-CPU',
	:p=>'skylake-himem',
        :autorun=>OPTIONS[:autorun],
        :excl=>" "}

## baits to do
def todo(v,i) 
      file = Dir.glob( i==1 ? "output/#{v}/totest.txt" : "output/#{v}/totest2.txt" )[0] 
      str=File.readlines(file)[0]
	if str == nil then
	return [];
end	
str.split(' ')
end

## baits done - rep i
def done(v,i)
  Dir.glob( "output/#{v}/rep-#{i}/*" ).map { |f|
    File.basename(f).sub!(".csv","")
  }
end

def brun(baits,v,i,args)
  if baits.length > 0 then
  fn = "#{v}-rep#{i}.txt"
  File.open(fn, "w") do |f|
    f.puts(baits)
  end
  n = (baits.length / 5).ceil
  # n = 20
  # args[:array] = "1-#{n}"
  # q=Qsub.new("#{v}-#{i}-#{COMFILE}", args)
  # q.add( "./do-tests.R --args patt=#{v} rep=#{i} ifile=#{fn} mult=10" )
  flag = OPTIONS[:autorun] ? "-r" : ""
  comm = "qR.rb #{flag} -c 1 -t \"8:00:00\" -j #{v}  -p skylake-himem -a CWALLACE-SL2-CPU -y 1-#{n} ./do-tests.R --args patt=#{v} rep=#{i} ifile=#{fn} mult=5"
  # comm = "qR.rb #{flag} -c 1 -t \"8:00:00\" -j #{v}  -a MRC-BSU-SL2-GPU -y 1-#{n} ./do-tests.R --args patt=#{v} rep=#{i} ifile=#{fn} mult=10"
  # system( comm )
  # comm = "qR.rb #{flag} -c 1 -t \"8:00:00\" -j #{v}#{i}  -y 1-#{n} ./do-tests.R --args patt=#{v} rep=#{i} ifile=#{fn} mult=10"
  system( comm )
  # system("qR.rb -c 1 -t \"6:00:00\" -j #{v} -y 1-#{n} ./do-tests.R --args patt=#{v} rep=#{i} ifile=#{fn} mult=10" )
  # q.close()
  end
end
  



################################################################################

## what cells to do
pchic = %w[ nCD4pchic aCD4pchic ]
val = %w[ aCD4val nCD4val ]
cells = []
if OPTIONS[:val]
  cells = cells + val
end
if OPTIONS[:pchic]
  cells = cells + pchic
end

puts cells
# cells= %w[ aCD4pchic ]

################################################################################

## define commands

if COMMAND=="summ" || COMMAND=="run" || COMMAND=="brun"
  cells.each { |v|
    # v="nCD4val"
    baits_todo=todo(v,1)
    baits_done1=done(v,1)
    baits_done2=done(v,2)
    baits_int = baits_done1 & baits_done2
    
    puts "#{v}, wanted : #{baits_todo.length()}"
    puts "\trep 1/2: #{baits_done1.length()} / #{baits_done2.length()} -> int #{baits_int.length()}"
    
    ## any 3/4 reps todo
    baits_todo2=todo(v,2)
    baits_done3=done(v,3)
    baits_done4=done(v,4)
    baits_int = baits_done3 & baits_done4
    
    puts "#{v}, wanted : #{baits_todo2.length()}"
    puts "\trep 3/4: #{baits_done3.length()} / #{baits_done4.length()} -> int #{baits_int.length()}"
    
    ## add to queue
    if COMMAND=="run"
      commands = commands + (baits_todo - baits_done1).map{ |b| "./do-tests.R --args patt=#{v} rep=1 bait=#{b}" }
      commands = commands + (baits_todo - baits_done2).map{ |b| "./do-tests.R --args patt=#{v} rep=2 bait=#{b}" }
      commands = commands + (baits_todo2 - baits_done3).map{ |b| "./do-tests.R --args patt=#{v} rep=3 bait=#{b}" }
      commands = commands + (baits_todo2 - baits_done4).map{ |b| "./do-tests.R --args patt=#{v} rep=4 bait=#{b}" }
    end
    
    if COMMAND=="brun"
      brun(baits_todo - baits_done1, v, 1, args) if baits_todo.length() > 0
      brun(baits_todo - baits_done2, v, 2, args) if baits_todo.length() > 0
      brun(baits_todo2 - baits_done3, v, 3, args) if baits_todo2.length() > 0
      brun(baits_todo2 - baits_done4, v, 4, args) if baits_todo2.length() > 0
    end
  }
end

if COMMAND == "corr" then
  args[:time] = "1:00:00"
  commands = cells.map { |v| "./corr.R --args patt=#{v}" }
end

if COMMAND == "mppc" then
  # args[:time] = "1:00:00"
  commands = cells.map { |v| "./mppc.R --args patt=#{v}" }
end 

if COMMAND == "clean" then
  system("rm *tmprun.sh*")
end

################################################################################

## run the things

if(commands.length > 0) then
  if OPTIONS[:int] then
    puts "RUNNING INTERACTIVELY"
    commands.each { |s|
      puts s
      system(s)
    }
  else
    nodes = (commands.length / (args[:tasks].to_f)).ceil
    ncomm = commands.length
    puts "  RUNNING #{ncomm} commands over #{nodes} nodes"
    q=Qsub.new(COMFILE, args)
    commands.each { |s| q.add(s) }
    q.close()
  end
else
  puts "  NO COMMANDS TO RUN"
end

